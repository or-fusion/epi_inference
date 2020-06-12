import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
import pandas as pd
from datetime import datetime
from . import reconstruction as recon
import pyutilib.misc.timing as timing
import json
import matplotlib.pyplot as plt

def run_multinode_mobility_time_varying_decay_lsq(cm_rep_cases, populations, mobility, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.

    Parameters
    ----------

    cm_rep_cases : a dataframe of *new* cases reported in
       each time period; each column in the dataframe is a separate time
       series
    populations : a dataframe with a single column that represents the
       population for different columns in cm_rep_cases
    sigma : float
       the rate constant for cases leaving the E compartment (1/incubation period)
    gamma : float
       the rate constant for leaving the I compartment. This will be multiplied
       by 3 when applied to each of the three I compartments
    report_delay : int
       the number of days between when someone is infected and when
       they will become a reported case (This should only shift the data
       and not impact the inference results.)
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    """
    # check the inputs
    assert sigma > 0
    assert gamma > 0
    assert report_delay > 0
    assert (populations > 0).all().all() == True
    assert reporting_factor >= 1

    for i in [-1]: #range(-6,2):
        # create the Pyomo optimization formulation
        regu = 1*10**i
        m = create_inference_regu_formulation(
            Cdates=Cdates,
            cumulative_reported_cases=cm_rep_cases,
            populations=populations,
            mobility=mobility,
            sigma=sigma,
            gamma_1=gamma*3,
            gamma_2=gamma*3,
            gamma_3=gamma*3,
            report_delay=report_delay,
            reporting_factor=reporting_factor,
            delta_beta_regu=regu,
            analysis_window=analysis_window,
            verbose=verbose
        )

        # call the solver
        solver = SolverFactory('ipopt')
        solver.options['tol']=1e-8
        status = solver.solve(m, tee=True) #verbose)
        m.display()

        # Check that the solve completed successfully
        if check_optimal_termination(status) == False:
            return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error.'}

        results = dict()
        results['date'] = [m.DATES[t] for t in m.TIMES]
        for i in m.NODES:
            betas = list()
            transmissions = 0
            for t in m.TIMES:
                transmissions += m.T_data[i][t]
                if transmissions < 20:
                    betas.append(None)
                else:
                    betas.append(pe.value(m.beta[i,t]))
            results[i] = betas

        df = pd.DataFrame(results)
        print(df)
        pd.set_option('display.max_rows', None)
        print(df.std())
        df.plot(x='date')
        plt.show()

        return {'results': results}

    """ OLD RESULTS STRUCTURE
    # build the dictionary of results
    betas = list()
    wdates = list()
    status = list()
    fips = list()
    pops = list()
    est_transmissions = list()
    window_days = list()
    for i in m.NODES:
        for w in m.WINDOWS:
            wdates.append(m.DATES[w])
            fips.append(i)
            pops.append(populations[i])
            est_transmissions.append(m.window_transmissions[i][w])
            window_days.append(m.window_days)
            if m.beta[i,w].stale == True or m.window_transmissions[i][w] <= 2:
                status.append('invalid_insufficient_data')
                betas.append(None)
            else:
                status.append('ok')
                betas.append(value(m.beta[i,w]))   

    ret =  {'dates': wdates, 'est_beta':betas, 'status':status, 'population': pops, 'est_transmissions':est_transmissions, 'window_days': window_days, 'FIPS':fips}
    #df = pd.DataFrame(ret)
    #df.to_csv('foo.csv')
    return ret
    """

def create_inference_regu_formulation(Cdates, cumulative_reported_cases, populations, mobility, sigma, gamma_1, gamma_2, gamma_3, report_delay, reporting_factor, delta_beta_regu, analysis_window, verbose=False):
    """
    Creates a one-step-ahead inference model using a decay
    model with 3 I compartments. The model is written in terms of absolute
    numbers of cases (not ln-transform).  The model combines estimates across
    multiple time series, one for each node.

    Parameters
    ----------
    Cdates: list of datetime objects
       The list of datetime objects that correspond to the dates for the
       cumulative_reported_cases
    cumulative_reported_cases : a dataframe of *new* cases reported in
       each time period; each column in the dataframe is a separate time
       series
    populations : a dataframe with a single column that represents the
       population for different columns in cumulative_reported_cases
    sigma : float
       the rate constant for cases leaving the E compartment (1/incubation period)
    gamma_1 : float
       the rate constant for leaving the I1 compartment.
    gamma_2 : float
       the rate constant for leaving the I2 compartment.
    gamma_3 : float
       the rate constant for leaving the I3 compartment.
    report_delay : int
       the number of days between when someone is infected and when
       they will become a reported case (This should only shift the data
       and not impact the inference results.)
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported

    """
    if len(analysis_window) != 0:
        raise NotImplementedError('analysis_window is not yet implemented for multinode_decay_lsq')
    model = pe.ConcreteModel()

    # Cached data
    model.sigma = sigma
    model.gamma_1 = gamma_1
    model.gamma_2 = gamma_2
    model.gamma_3 = gamma_3
    model.eta = 0.5 # fraction of the day spent "away"
    model.report_delay = report_delay
    model.reporting_factor = reporting_factor
    model.delta_beta_regu = delta_beta_regu

    #model.NODES = pe.Set(initialize=list(range(len(cumulative_reported_cases.keys()))))
    model.NODES = pe.Set(initialize=list(k for k in cumulative_reported_cases.keys()))

    model.mobility = dict(mobility)
    model.MOBILITY_TUPLES = list()
    #model.mobility = dict()
    for i in model.NODES:
        if i not in model.mobility:
            model.mobility[i] = {}
        for j in model.mobility[i]:
            if i in model.NODES and j in model.NODES:
                model.MOBILITY_TUPLES.append((i,j))
    model.populations = dict()
    
    model.T_data = dict()
    model.I_data = dict()
    model.S_data = dict()
    for nodeid in model.NODES:
        model.populations[nodeid] = float(populations[nodeid]) # overall population
        cm_rep_cases_node = [v for v in cumulative_reported_cases[nodeid].values]

        rdates, rcases, dates, T, S, E, I1, I2, I3, R = \
            recon.reconstruct_states_deterministic_decay(
                Cdates=Cdates,
                cumulative_reported_cases=cm_rep_cases_node,
                population=model.populations[nodeid],
                sigma=sigma,
                gamma=gamma_1/3,
                reporting_factor=reporting_factor,
                report_delay=report_delay
            )
        model.T_data[nodeid] = T
        model.I_data[nodeid] = dict()
        model.I_data[nodeid]['I1'] = I1
        model.I_data[nodeid]['I2'] = I2
        model.I_data[nodeid]['I3'] = I3
        model.S_data[nodeid] = S
        
        if not hasattr(model, 'TIMES'):
            model.TIMES = pe.Set(initialize=[i for i in range(len(T))], ordered=True)
        if not hasattr(model, 'DATES'):
            model.DATES = dates

    model.beta = pe.Var(model.NODES, model.TIMES, initialize=1.0, bounds=(0,None)) # transmission parameter
    # for now, alpha is not used
    # model.alpha = pe.Var(initialize=1.0)
    # model.alpha.fix(1.0)

    # define the variable for estimated transmissions
    model.T_hat = pe.Var(model.NODES, model.TIMES, initialize=1.0)
    # infection process
    def _infection_process(m, i, t):
        percent_mobile = 0
        if i in m.mobility:
            percent_mobile = sum(m.mobility[i][j]/m.populations[i] for j in m.mobility[i] if j in m.NODES)

        return m.T_hat[i,t] == m.beta[i,t] * (m.I_data[i]['I1'][t] + m.I_data[i]['I2'][t] + m.I_data[i]['I3'][t]) / m.populations[i] * m.S_data[i][t] * (1-m.eta*percent_mobile) \
            + sum(m.beta[j,t] * (m.I_data[j]['I1'][t] + m.I_data[j]['I2'][t] + m.I_data[j]['I3'][t]) / m.populations[j] * m.S_data[i][t] * (m.eta*m.mobility[i][j]/m.populations[i]) for j in m.mobility[i] if j in m.NODES)

    model.infection_process = pe.Constraint(model.NODES, model.TIMES, rule=_infection_process)

    model.delta_beta = pe.Var(model.NODES, model.TIMES, initialize=0)
    def _delta_beta_con(m, i, t):
        if t == m.TIMES.first():
            return pe.Constraint.Skip
        return m.delta_beta[i,t] == m.beta[i,t] - m.beta[i,t-1]
    model.delta_beta_con = pe.Constraint(model.NODES, model.TIMES, rule=_delta_beta_con)

    # least squares objective function
    def _lse(m):
        return sum( (m.T_hat[i,t] - m.T_data[i][t])**2 for i in m.NODES for t in m.TIMES)
    model.lse = pe.Expression(rule=_lse)

    def _regu(m):
        return sum(m.delta_beta[i,t]**2 for i in m.NODES for t in m.TIMES)
    model.regu = pe.Expression(rule=_regu)
    
    def _total_lse(m):
        return m.lse + m.delta_beta_regu * m.regu
    model.total_lse = pe.Objective(rule=_total_lse)

    return model

