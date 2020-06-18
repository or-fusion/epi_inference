__all__ = ['run_multinode_mobility_window_decay_lsq_iterative']

import pyutilib.misc.timing as timing
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination

from .util import get_windows

def run_multinode_mobility_window_decay_lsq_iterative(*, recon, mobility, analysis_window, select_window=None, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.

    Parameters
    ----------
    recon : dict()
       A dictionary with reconstruction data, indexed by FIPS codes for US counties.
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    """
    timing.tic('Starting timer for model construction - Pyomo')
    #
    # Create the Pyomo optimization formulation
    #
    window = int(analysis_window.get('days',14))
    assert(window >= 1)

    eta = 0.5 # fraction of the day spent "away"
    nodes = set(k for k in recon)

    T_data = dict()
    I1_data = dict()
    I2_data = dict()
    I3_data = dict()
    S_data = dict()
    populations = dict()
    percent_mobile = dict()
    dates = None
    for nodeid in nodes:
        T_data[nodeid] = recon[nodeid]['transmissions']
        I1_data[nodeid] = recon[nodeid]['I1']
        I2_data[nodeid] = recon[nodeid]['I2']
        I3_data[nodeid] = recon[nodeid]['I3']
        S_data[nodeid] = recon[nodeid]['S']
        populations[nodeid] = recon[nodeid]['population']
        percent_mobile[nodeid] = sum(mobility[nodeid][j] for j in mobility[nodeid] if j in nodes)/populations[nodeid] if nodeid in mobility else 0

        if dates is None:
            dates = recon[nodeid]['dates']
    timing.toc('setup population and mobility information')

    # define the tuples for the windows
    windows = get_windows(dates, window_days=window, select_window=select_window)
    WINDOW_TIMES = windows.WINDOW_TIMES

    # get the approximate transmissions over the window period
    window_transmissions = dict()
    for i in nodes:
        d = dict()
        for w in WINDOW_TIMES:
            d[w] = sum(T_data[i][t] for t in WINDOW_TIMES[w])
        window_transmissions[i] = d
        
    # Setup results object
    results = {}
    for i in recon:
        county = {}
        county['FIPS'] = i
        county['window_days'] = window
        county['date'] = [recon[i]['dates'][w] for w in WINDOW_TIMES]
        if i in nodes:
            county['population'] = recon[i]['population']
            county['beta'] = []
            county['status'] = []
            county['infections_in_window'] = []
        results[i] = county

    #
    # Setup and solve different problems for each window
    #
    for w in WINDOW_TIMES:
        timing.tic('Starting timer for model construction - Pyomo')
        model = pe.ConcreteModel()
        model.NODES = pe.Set(initialize=nodes, ordered=False)
        model.beta = pe.Var(model.NODES, initialize=1.0, bounds=(0,None)) 

        infections = 0
        for t in WINDOW_TIMES[w]:
            for i in nodes:
                infections += I1_data[i][t] + I2_data[i][t] + I3_data[i][t]
        if infections > 0:
            # define the expression for estimated transmissions
            def _infection_process(m, i, t):
                return model.beta[i] * ((I1_data[i][t] + I2_data[i][t] + I3_data[i][t]) /populations[i] * S_data[i][t] * (1-eta*percent_mobile[i])) + sum(model.beta[j] * ((I1_data[j][t] + I2_data[j][t] + I3_data[j][t]) * S_data[i][t] * mobility[i][j] * eta / (populations[j]*populations[i])) for j in mobility[i] if j in nodes)
    
            model.T_hat = pe.Expression(model.NODES, WINDOW_TIMES[w], rule=_infection_process)
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

        # Grab the results from the solver
        for i in recon:
            if i not in nodes:
                continue
            county = results[i]

            if model.beta[i].stale == True:
                county['beta'].append( None )
                county['status'].append( 'stale' )
            else:
                county['beta'].append( value(model.beta[i]) )
                county['status'].append( 'ok' )
            county['infections_in_window'].append( window_transmissions[i][w] )

    return results

