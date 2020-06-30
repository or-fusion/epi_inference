__all__ = ['run_multinode_mobility_window_decay_lsq_iterative']

import pyutilib.misc.timing as timing
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination

from .util import get_windows, indices_since_first_nonzero

def run_multinode_mobility_window_decay_lsq_iterative(*, recon, mobility, analysis_window, objective, select_window=None, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.

    Parameters
    ----------
    recon : dict()
       A dictionary with reconstruction data, indexed by FIPS codes for US counties.
    mobility: dict()
       A dictionary of sparse mobility data
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the inference. If None, then the full set of data will be used.
       The key "days" indicates the length of the analysis window in days. The
       date returned in the analysis is the last day of this window.
    objective : str
       Choice of objective function to use. Least-squares (lsq) is the default.
    select_window: str
       Date corresponding to a particular window for analysis
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    """
    timing.tic('Starting timer for model construction - Pyomo')
    # check the input data
    assert objective == 'lsq' or objective == 'binomial-likelihood'
    window_length_days = int(analysis_window.get('days',14))
    assert(window_length_days >= 1)

    eta = 0.5 # fraction of the day spent "away"

    # create the set of nodes used in the model - use sorted_nodes for
    # determinism, the set nodes for check
    nodes = set(k for k in recon)
    sorted_nodes = sorted(nodes)

    # build data structures for the parameters
    T_data = dict()
    I1_data = dict()
    I2_data = dict()
    I3_data = dict()
    S_data = dict()
    orig_rep_cases = dict()
    days_since_first_reported = dict()
    populations = dict()
    percent_mobile = dict()
    dates = None
    for nodeid in sorted_nodes:
        T_data[nodeid] = recon[nodeid]['transmissions']
        I1_data[nodeid] = recon[nodeid]['I1']
        I2_data[nodeid] = recon[nodeid]['I2']
        I3_data[nodeid] = recon[nodeid]['I3']
        S_data[nodeid] = recon[nodeid]['S']
        orig_rep_cases[nodeid] = recon[nodeid]['orig_rep_cases']
        days_since_first_reported[nodeid] = indices_since_first_nonzero(orig_rep_cases[nodeid])
        populations[nodeid] = recon[nodeid]['population']
        percent_mobile[nodeid] = sum(mobility[nodeid][j] for j in mobility[nodeid] if j in nodes)/populations[nodeid] if nodeid in mobility else 0

        if dates is None:
            # dates should be the same across all counties
            dates = recon[nodeid]['dates']
        assert dates == recon[nodeid]['dates']

    timing.toc('setup inference parameters (transmissions, population, mobility, etc.')

    # define the WINDOW_TIMES tuple pairs for the windows
    windows = get_windows(dates, window_days=window_length_days, select_window=select_window)
    WINDOW_TIMES = windows.WINDOW_TIMES

    # gather some extra data to be reported with the window period
    # approx. transmissions and infections (from reconstruction)
    # and reported cases,
    window_transmissions = dict()
    infectious_pop_over_window = dict()
    transmissions_over_window = dict()
    reported_cases_over_window = dict()
    days_since_first_reported_by_window = dict()
    for i in sorted_nodes:
        infectious_pop_over_window[i] = dict()
        transmissions_over_window[i] = dict()
        reported_cases_over_window[i] = dict()
        days_since_first_reported_by_window[i] = dict()
        d = dict()
        for w in WINDOW_TIMES:
            d[w] = sum(T_data[i][t] for t in WINDOW_TIMES[w])
            infectious_pop_over_window[i][w] = list(I1_data[i][t] + I2_data[i][t] + I3_data[i][t] for t in WINDOW_TIMES[w])
            transmissions_over_window[i][w] = list(T_data[i][t] for t in WINDOW_TIMES[w])
            reported_cases_over_window[i][w] = list(orig_rep_cases[i][t] for t in WINDOW_TIMES[w])
            days_since_first_reported_by_window[i][w] = days_since_first_reported[i][WINDOW_TIMES[w][-1]]
        window_transmissions[i] = d

    # Setup results object
    results = {}
    for i in recon:
        county = {}
        county['FIPS'] = i
        county['window_days'] = window_length_days
        county['date'] = [recon[i]['dates'][w] for w in WINDOW_TIMES]
        if i in nodes:
            county['population'] = recon[i]['population']
            county['beta'] = []
            county['status'] = []
            county['infections_in_window'] = []
            county['infectious_pop_over_window'] = []
            county['transmissions_over_window'] = []
            county['reported_cases_over_window'] = []
            county['days_since_first_reported'] = []
        results[i] = county

    #
    # Setup and solve different problems for each window
    #
    for w in WINDOW_TIMES:
        timing.tic('Starting timer for model construction - Pyomo')
        model = pe.ConcreteModel()
        model.NODES = pe.Set(initialize=sorted_nodes, ordered=True)
        model.beta = pe.Var(model.NODES, initialize=1.0, bounds=(0,None)) 

        # check the total number of infections - if there are none across
        # all counties, the optimization will not run
        total_infections = 0
        for t in WINDOW_TIMES[w]:
            for i in sorted_nodes:
                total_infections += I1_data[i][t] + I2_data[i][t] + I3_data[i][t]

        if total_infections > 0:
            # define the expression for estimated transmissions
            def _infection_process(m, i, t):
                return model.beta[i] * ((I1_data[i][t] + I2_data[i][t] + I3_data[i][t]) /populations[i] * S_data[i][t] * (1-eta*percent_mobile[i])) + sum(model.beta[j] * ((I1_data[j][t] + I2_data[j][t] + I3_data[j][t]) * S_data[i][t] * mobility[i][j] * eta / (populations[j]*populations[i])) for j in mobility[i] if j in nodes)
    
            model.T_hat = pe.Expression(model.NODES, WINDOW_TIMES[w], rule=_infection_process)

            # we want to record which T_hat expressions are guaranteed to be zero
            # so we can remove them from the likelihood function (don't want ln(0))
            # since they are initialized to 1.0, we can easily check this by
            # evaluating them. If they evaluate to zero, then we record the indices
            zero_T_hats = set()
            for n in model.NODES:
                for t in WINDOW_TIMES[w]:
                    if abs(pe.value(model.T_hat[n,t])) < 1e-10:
                        zero_T_hats.add((n,t))

            timing.toc('built infection process')

            model.total_lse = pe.Objective(expr=sum((model.T_hat[i,t] - T_data[i][t])**2 for i,t in model.T_hat)/len(model.T_hat))
            timing.toc('built objective')

            # call the solver
            timing.tic('Starting timer for solver')
            solver = SolverFactory('ipopt')
            solver.options['tol']=1e-8
            status = solver.solve(model, options={'print_level':0})
            timing.toc('Finished solver')

            # Check that the solve completed successfully
            if check_optimal_termination(status) == False:
                return {'beta': None, 'status': 'failed', 'msg': 'Unknown solver error for window %s and time %s.' % (str(w),str(t))}

            if objective == 'binomial-likelihood':
                # solve the binomial-likelihood
                # note that we want the call to the lsq part above to initialize
                model.total_like = pe.Objective(
                    expr=sum( T_data[i][t]*pe.log(1.0-pe.exp(-model.T_hat[i,t]/S_data[i][t])) + (S_data[i][t]-T_data[i][t])*(-model.T_hat[i,t]/S_data[i][t]) \
                              for i,t in model.T_hat if (i,t) not in zero_T_hats),
                    sense=pe.maximize
                )
                model.total_lse.deactivate()
                timing.tic('Starting timer for solver with likelihood objective')
                status = solver.solve(model, tee=True)
                timing.toc('Finished solver with likelihood objective')

                # Check that the solve completed successfully
                if check_optimal_termination(status) == False:
                    return {'beta': None, 'status': 'failed-likelihood', 'msg': 'Unknown solver error for window %s and time %s.' % (str(w),str(t))}

        # Grab the results from the solver
        for i in recon:
            if i not in nodes:
                continue
            county = results[i]

            if total_infections <= 0 or model.beta[i].stale == True:
                # did not have sufficient data to determine a value for beta
                county['beta'].append( None )
                county['status'].append( 'stale' )
            else:
                county['beta'].append( value(model.beta[i]) )
                county['status'].append( 'ok' )
            county['infections_in_window'].append( window_transmissions[i][w] )
            county['infectious_pop_over_window'].append( infectious_pop_over_window[i][w] )
            county['transmissions_over_window'].append( transmissions_over_window[i][w] )
            county['reported_cases_over_window'].append( reported_cases_over_window[i][w] )
            county['days_since_first_reported'].append( days_since_first_reported_by_window[i][w] )

    return results
