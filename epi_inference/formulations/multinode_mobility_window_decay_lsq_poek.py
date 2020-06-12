__all__ = ['run_multinode_mobility_window_decay_lsq_poek']

import itertools
import pyutilib.misc.timing as timing
from pyutilib.misc.misc import Bunch
try:
    import poek as pk
    poek_available = True
except:
    poek_available = False


def run_multinode_mobility_window_decay_lsq_poek(*, recon, mobility, analysis_window, verbose=False):
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
    if not poek_available:
        raise RuntimeError("Cannot solve the mobility window formulation with poek")

    # create the Pyomo optimization formulation
    m = create_inference_window_formulation(
        recon=recon,
        mobility=mobility,
        analysis_window=analysis_window,
        verbose=verbose
    )

    if m.model is None:
        return {'beta': None, 'status': 'failed', 'msg': 'Empty model.'}

    #m.model.write("output.nl")
    # call the solver
    timing.tic('Starting timer for solver')
    nlp = pk.nlp_model(m.model, "cppad")
    solver = pk.nlp_solver('ipopt')
    #solver.options['tol']=1e-8
    status = solver.solve(nlp)
    timing.toc('Finished solver')

    # Check that the solve completed successfully
    ##if check_optimal_termination(status) == False:
    ##    return {'beta': None, 'status': 'failed', 'msg': 'Unknown solver error.'}

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
                if False and  m.beta[i,w].stale == True:
                    county['beta'].append( None )
                    county['status'].append( 'stale' )
                else:
                    county['beta'].append( m.beta[i,w].value )
                    county['status'].append( 'ok' )
                county['infections_in_window'].append( m.window_transmissions[i][w] )
        results[i] = county

    return results


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

    timing.tic('Starting timer for model construction - POEK')
    model = pk.model()

    eta = 0.5 # fraction of the day spent "away"

    NODES = set(k for k in recon)

    T_data = dict()
    I1_data = dict()
    I2_data = dict()
    I3_data = dict()
    S_data = dict()
    populations = dict()
    percent_mobile = dict()
    for nodeid in NODES:
        T_data[nodeid] = recon[nodeid]['transmissions']
        I1_data[nodeid] = recon[nodeid]['I1']
        I2_data[nodeid] = recon[nodeid]['I2']
        I3_data[nodeid] = recon[nodeid]['I3']
        S_data[nodeid] = recon[nodeid]['S']
        populations[nodeid] = recon[nodeid]['population']
        percent_mobile[nodeid] = sum(mobility[nodeid][j] for j in mobility[nodeid] if j in NODES)/populations[nodeid] if nodeid in mobility else 0

        if not hasattr(model, 'TIMES'):
            TIMES = [i for i in range(len(recon[nodeid]['transmissions']))]
    timing.toc('setup population and mobility information')

    # define the tuples for the windows
    WINDOWS = list()
    WINDOW_TIMES = list()
    window_days = window
    for i in range(len(TIMES)):
        if i % 7 != 0:
            continue
        if i < window_days:
            continue
        for j in range(i+1-window_days, i+1):
            WINDOW_TIMES.append((i,j)) 
        WINDOWS.append(i)
    timing.toc('built windows')

    # transmission parameter
    beta = model.variable(index=[(i,j) for i in NODES for j in WINDOWS], value=1.0, lb=0)
    timing.toc('built variables')

    T_hat = {}
    for i,W in itertools.product(NODES, WINDOW_TIMES):
            w,t = W
    #for i in NODES:
    #    for w,t in WINDOW_TIMES:

            T_hat[i,w,t] = beta[i,w] * ((I1_data[i][t] + I2_data[i][t] + I3_data[i][t]) /populations[i] * S_data[i][t] * (1-eta*percent_mobile[i]))  + sum(beta[j,w] * ((I1_data[j][t] + I2_data[j][t] + I3_data[j][t]) * S_data[i][t] * mobility[i][j] * eta / (populations[i]*populations[j])) for j in mobility[i] if j in NODES)

    timing.toc('built infection process')

    # least squares objective function
    lse = {}
    for i in NODES:
        lse[i] = sum( (T_hat[i,w,t] - T_data[i][t])**2 for w,t in WINDOW_TIMES)

    model.add( sum(lse[i] for i in NODES) )
    timing.toc('built objective')

    # get the approximate transmissions over the window period
    window_transmissions = dict()
    for i in NODES:
        d = dict()
        for w in WINDOWS:
            d[w] = sum(T_data[i][t] for ww,t in WINDOW_TIMES if ww == w)
        window_transmissions[i] = d
        
    return Bunch(model=model, WINDOWS=WINDOWS, window_transmissions=window_transmissions, beta=beta, NODES=NODES, window_days=window_days)

