__all__ = ['run_multinode_mobility_window_decay_lsq']

import pyutilib.misc.timing as timing
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination

from .util import get_windows


def run_multinode_mobility_window_decay_lsq(*, recon, mobility, analysis_window, select_window=None, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.

    Parameters
    ----------
    recon : dict()
       A dictionary with reconstruction data, indexed by FIPS codes for US counties.
    mobility : dict()
        A dictionary with inter-county mobility rates.
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    select_window : str
       ISO date format that the window that will be used in this estimation.  If None,
       then all windows are used.
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    """
    # create the Pyomo optimization formulation
    m = create_inference_window_formulation(
        recon=recon,
        mobility=mobility,
        analysis_window=analysis_window,
        select_window=select_window,
        verbose=verbose
    )

    if m is None:
        return {'beta': None, 'status': 'failed', 'msg': 'Empty model.'}

    # call the solver
    timing.tic('Starting timer for solver')
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m)
    timing.toc('Finished solver')

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'beta': None, 'status': 'failed', 'msg': 'Unknown solver error.'}

    results = {}
    for i in recon:
        county = {}
        county['FIPS'] = i
        county['window_days'] = m.window_days
        county['date'] = [recon[i]['dates'][w] for w in m.WINDOWS]
        if i in m.NODES:
            county['population'] = recon[i]['population']
            county['beta'] = []
            county['status'] = []
            county['infections_in_window'] = []
            for w in m.WINDOWS:
                if m.beta[i,w].stale == True:
                    county['beta'].append( None )
                    county['status'].append( 'stale' )
                else:
                    county['beta'].append( value(m.beta[i,w]) )
                    county['status'].append( 'ok' )
                county['infections_in_window'].append( m.window_transmissions[i][w] )
        results[i] = county

    return results


def create_inference_window_formulation(*, recon, mobility, analysis_window, select_window, verbose=False):
    """
    Creates a one-step-ahead inference model using a decay
    model with 3 I compartments. The model is written in terms of absolute
    numbers of cases (not ln-transform).  The model combines estimates across
    multiple time series, one for each node.

    Parameters
    ----------
    recon : dict()
       A dictionary with reconstruction data, indexed by FIPS codes for US counties.
    mobility : dict()
        A dictionary with inter-county mobility rates.
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    """
    window = int(analysis_window.get('days',14))
    assert(window >= 1)

    timing.tic('Starting timer for model construction - Pyomo')
    model = pe.ConcreteModel()

    eta = 0.5 # fraction of the day spent "away"

    nodes = set(k for k in recon)
    model.NODES = pe.Set(initialize=nodes)

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
    model.TIMES = pe.Set(initialize=windows.TIMES, ordered=True)
    WINDOWS = windows.WINDOWS
    WINDOW_TIMES = windows.WINDOW_TIMES_LIST

    # transmission parameter
    model.beta = pe.Var(model.NODES, WINDOWS, initialize=1.0, bounds=(0,None)) 
    timing.toc('built variables')

    # define the expression for estimated transmissions
    def _infection_process(m, i, w, t):

        return m.beta[i,w] * ((I1_data[i][t] + I2_data[i][t] + I3_data[i][t]) /populations[i] * S_data[i][t] * (1-eta*percent_mobile[i])) + sum(m.beta[j,w] * ((I1_data[j][t] + I2_data[j][t] + I3_data[j][t]) * S_data[i][t] * mobility[i][j] * eta / (populations[j]*populations[i])) for j in mobility[i] if j in nodes)

    model.T_hat = pe.Expression(model.NODES, WINDOW_TIMES, rule=_infection_process)
    timing.toc('built infection process')

    # least squares objective function
    def _lse(m, i):
        return sum( (m.T_hat[i,w,t] - T_data[i][t])**2 for w,t in WINDOW_TIMES) #\
    model.lse = pe.Expression(model.NODES, rule=_lse)

    model.total_lse = pe.Objective(expr=sum(model.lse[i] for i in nodes))
    timing.toc('built objective')

    # get the approximate transmissions over the window period
    model.window_transmissions = dict()
    for i in nodes:
        d = dict()
        for w in WINDOWS:
            d[w] = sum(T_data[i][t] for ww,t in WINDOW_TIMES if ww == w)
        model.window_transmissions[i] = d
        
    model.WINDOWS = WINDOWS
    model.window_days = window
    return model

