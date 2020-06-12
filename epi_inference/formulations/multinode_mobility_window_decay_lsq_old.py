__all__ = ['run_multinode_mobility_window_decay_lsq_old']

import pyutilib.misc.timing as timing
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination


def run_multinode_mobility_window_decay_lsq_old(*, recon, mobility, analysis_window, verbose=False):
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
    # create the Pyomo optimization formulation
    m = create_inference_window_formulation(
        recon=recon,
        mobility=mobility,
        analysis_window=analysis_window,
        verbose=verbose
    )

    # call the solver
    timing.tic('Starting timer for solver')
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m)
    timing.toc('Finished solver')

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error.'}

    results = list()
    for i in m.NODES:
        for w in m.WINDOWS:
            d = dict()
            d['date'] = m.DATES[w]
            d['window_days'] = m.window_days
            if m.beta[i,w].stale == True:
                d['est_beta'] = None
                d['status'] = 'stale'
            else:
                d['est_beta'] = value(m.beta[i,w])
                d['status'] = 'ok'
            d['population'] = recon[i]['population']
            d['infections_in_window'] = m.window_transmissions[i][w]
            d['FIPS'] = i
            results.append(d)

    return {'results': results}


def create_inference_window_formulation(*, recon, mobility, analysis_window, verbose=False):
    """
    Creates a one-step-ahead inference model using a decay
    model with 3 I compartments. The model is written in terms of absolute
    numbers of cases (not ln-transform).  The model combines estimates across
    multiple time series, one for each node.

    Parameters
    ----------
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.

    """
    window = int(analysis_window.get('days',14))
    assert(window >= 1)

    timing.tic('Starting timer for model construction')
    model = pe.ConcreteModel()

    # Cached data
    #model.sigma = sigma
    #model.gamma_1 = gamma_1
    #model.gamma_2 = gamma_2
    #model.gamma_3 = gamma_3
    model.eta = 0.5 # fraction of the day spent "away"
    #model.report_delay = report_delay
    #model.reporting_factor = reporting_factor
    #model.delta_beta_regu = 1e-4

    model.NODES = pe.Set(initialize=list(k for k in recon.keys()))

    model.mobility = dict(mobility)
    model.MOBILITY_TUPLES = list()
    for i in model.NODES:
        if i not in model.mobility:
            model.mobility[i] = {}
        for j in model.mobility[i]:
            if i in model.NODES and j in model.NODES:
                model.MOBILITY_TUPLES.append((i,j))
    
    model.T_data = dict()
    model.I_data = dict()
    model.S_data = dict()
    model.populations = dict()
    for nodeid in model.NODES:
        #model.populations[nodeid] = float(populations[nodeid]) # overall population
        #cm_rep_cases_node = [v for v in cumulative_reported_cases[nodeid].values]
        #
        #rdates, rcases, dates, T, S, E, I1, I2, I3, R = \
        #    recon.reconstruct_states_deterministic_decay(
        #        Cdates=Cdates,
        #        cumulative_reported_cases=cm_rep_cases_node,
        #        population=model.populations[nodeid],
        #        sigma=sigma,
        #        gamma=gamma_1/3,
        #        reporting_factor=reporting_factor,
        #        report_delay=report_delay
        #    )
        model.T_data[nodeid] = recon[nodeid]['transmissions']
        model.I_data[nodeid] = dict()
        model.I_data[nodeid]['I1'] = recon[nodeid]['I1']
        model.I_data[nodeid]['I2'] = recon[nodeid]['I2']
        model.I_data[nodeid]['I3'] = recon[nodeid]['I3']
        model.S_data[nodeid] = recon[nodeid]['S']
        
        model.populations[nodeid] = recon[nodeid]['population']

        if not hasattr(model, 'TIMES'):
            model.TIMES = pe.Set(initialize=[i for i in range(len(recon[nodeid]['transmissions']))], ordered=True)
        if not hasattr(model, 'DATES'):
            model.DATES = recon[nodeid]["dates"]
    timing.toc('setup population and mobility information')

    # define the tuples for the windows
    model.WINDOWS = list()
    model.WINDOW_TIMES = list()
    model.window_days = window
    for i in range(len(model.TIMES)):
        if i % 7 != 0:
            continue
        if i < model.window_days:
            continue
        for j in range(i+1-model.window_days, i+1):
            model.WINDOW_TIMES.append((i,j)) 
        model.WINDOWS.append(i)
    timing.toc('built windows')

    # get the approximate transmissions over the window period
    model.window_transmissions = dict()
    for i in model.NODES:
        d = dict()
        for w in model.WINDOWS:
            d[w] = sum(model.T_data[i][t] for ww,t in model.WINDOW_TIMES if ww == w)
        model.window_transmissions[i] = d
    #print(WINDOWS)
    #for t in WINDOW_TIMES:
    #    print(t)
    #quit()
        
    model.beta = pe.Var(model.NODES, model.WINDOWS, initialize=1.0, bounds=(0,None)) # transmission parameter
    # for now, alpha is not used
    # model.alpha = pe.Var(initialize=1.0)
    # model.alpha.fix(1.0)

    # define the variable for estimated transmissions
    model.T_hat = pe.Var(model.NODES, model.WINDOW_TIMES, initialize=1.0)
    timing.toc('built variables')
    # infection process
    def _infection_process(m, i, w, t):
        percent_mobile = 0
        if i in m.mobility:
            percent_mobile = sum(m.mobility[i][j]/m.populations[i] for j in m.mobility[i] if j in m.NODES)

        return m.T_hat[i,w,t] == m.beta[i,w] * (m.I_data[i]['I1'][t] + m.I_data[i]['I2'][t] + m.I_data[i]['I3'][t]) / m.populations[i] * m.S_data[i][t] * (1-m.eta*percent_mobile) \
            + sum(m.beta[j,w] * (m.I_data[j]['I1'][t] + m.I_data[j]['I2'][t] + m.I_data[j]['I3'][t]) / m.populations[j] * m.S_data[i][t] * (m.eta*m.mobility[i][j]/m.populations[i]) for j in m.mobility[i] if j in m.NODES)

    model.infection_process = pe.Constraint(model.NODES, model.WINDOW_TIMES, rule=_infection_process)
    timing.toc('built infection process')
    """
    model.delta_beta = pe.Var(model.MOBILITY_TUPLES, model.WINDOWS)
    def _delta_beta_con(m, i, j, w):
        return m.delta_beta[i,j,w] == m.beta[i,w] - m.beta[j,w]
    model.delta_beta_con = pe.Constraint(model.MOBILITY_TUPLES, model.WINDOWS, rule=_delta_beta_con)
    """
    # least squares objective function
    def _lse(m, i):
        return sum( (m.T_hat[i,w,t] - m.T_data[i][t])**2 for w,t in m.WINDOW_TIMES) #\
            #+ m.delta_beta_regu*sum( m.mobility[i][j]*m.delta_beta[i,j,w]**2 for j in m.mobility[i] for w in m.WINDOWS if i in m.NODES and j in m.NODES)
    model.lse = pe.Expression(model.NODES, rule=_lse)

    def _total_lse(m):
        return sum( m.lse[i] for i in m.NODES )
    model.total_lse = pe.Objective(rule=_total_lse)
    timing.toc('built objective')

    return model

