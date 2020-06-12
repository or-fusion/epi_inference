import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
import pandas as pd
from datetime import datetime
import epi_inference.reconstruction.common as rcommon
import epi_inference.reconstruction.deterministic as recond

def run_multinode_mobility_decay_lsq(cm_rep_cases, populations, mobility, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
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

    # ToDo: this needs to be passed in - for now create a default set of dates
    #Cdates = pd.date_range(end=datetime(year=2020, month=4, day=12), periods=len(cm_rep_cases)).to_pydatetime().tolist()

    # create the Pyomo optimization formulation
    m = create_inference_formulation_multinode_mobility_decay_lsq(
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
        analysis_window=analysis_window,
        verbose=verbose
    )

    # call the solver
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error.'}

    results = list()
    for i in m.NODES:
        d = dict()
        if m.beta[i].stale == True:
            d['est_beta'] = None
            d['status'] = 'stale'
        else:
            d['est_beta'] = value(m.beta[i])
            d['status'] = 'ok'
        d['population'] = populations[i]
        d['total_cases'] = cm_rep_cases[i][-1]
        d['FIPS'] = i
        results.append(d)

    return {'results': results}

    """ OLD RESULTS STRUCTURE
    # check that the beta value was successfully solved
    betas = list()
    status = list()
    fips = list()
    pops = list()
    casecount = list()
    for i in m.NODES:
        fips.append(i)
        pops.append(populations[i])
        casecount.append(cm_rep_cases[i][-1])
        if m.beta[i].stale == True or cm_rep_cases[i][-1] < 3:
            status.append('invalid_insufficient_data')
            betas.append(None)
        else:
            status.append('ok')
            betas.append(value(m.beta[i]))   

    return {'est_beta':betas, 'status':status, 'population': pops, 'case_count':casecount, 'FIPS':fips}
    """

    """
    for i in m.beta:
        print('{},{},{},{},{}'.format(i, value(m.beta[i]), float(cm_rep_cases[i].values[-1]), populations[i], m.beta[i].stale))

    import matplotlib.pyplot as plt
    df = pd.DataFrame({'est_beta': [value(m.beta[i]) for i in m.beta if float(cm_rep_cases[i].values[-1]) > 10]})
    df.hist()
    plt.show()

    errors = dict()
    for c in m.NODES:
        cerrors = list()
        for t in m.TIMES:
            cerrors.append( (value(m.T_hat[c,t])-m.T_data[c][t])/max(1.0,m.T_data[c][t]) )
        errors[c] = cerrors

    df = pd.DataFrame(errors)
    df.plot()
    plt.show()
    quit()


    return {'est_beta': value(m.beta),'status': 'ok', 'msg': 'Optimal solution found'}
    """


def create_inference_formulation_multinode_mobility_decay_lsq(Cdates, cumulative_reported_cases, populations, mobility, sigma, gamma_1, gamma_2, gamma_3, report_delay, reporting_factor, analysis_window, verbose=False):
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

    #model.NODES = pe.Set(initialize=list(range(len(cumulative_reported_cases.keys()))))
    model.NODES = pe.Set(initialize=list(k for k in cumulative_reported_cases.keys()))
    model.beta = pe.Var(model.NODES, initialize=1.0, bounds=(0,None)) # transmission parameter
    # for now, alpha is not used
    # model.alpha = pe.Var(initialize=1.0)
    # model.alpha.fix(1.0)

    model.mobility = dict(mobility)
    #model.mobility = dict()
    for i in model.NODES:
        if i not in model.mobility:
            model.mobility[i] = {}
    model.populations = dict()
    
    model.T_data = dict()
    model.I_data = dict()
    model.S_data = dict()
    for nodeid in model.NODES:
        model.populations[nodeid] = float(populations[nodeid]) # overall population
        cm_rep_cases_node = [v for v in cumulative_reported_cases[nodeid].values]
        rdates, rcases = \
            rcommon.reported_cases_from_cumulative(Cdates, cm_rep_cases_node)

        dates, T, S, E, I1, I2, I3, R = \
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

    # define the variable for estimated transmissions
    model.T_hat = pe.Var(model.NODES, model.TIMES, initialize=1.0)

    # infection process
    #
    # NOTE: We need to be careful how we model borders.  If the mobility matrix says we cross them, then 
    #   how do we model the impact of infection on that mobile population?
    #
    def _infection_process(m, i, t):
        percent_mobile = 0
        if i in m.mobility:
            percent_mobile = sum(m.mobility[i][j]/m.populations[i] for j in m.mobility[i] if j in m.NODES)
            #percent_mobile = sum(m.mobility[i][j]/m.populations[i] for j in m.mobility[i])
        return m.T_hat[i,t] == m.beta[i] * (m.I_data[i]['I1'][t] + m.I_data[i]['I2'][t] + m.I_data[i]['I3'][t]) / m.populations[i] * m.S_data[i][t] * (1-m.eta*percent_mobile) \
            + sum(m.beta[j] * (m.I_data[j]['I1'][t] + m.I_data[j]['I2'][t] + m.I_data[j]['I3'][t]) / m.populations[j] * m.S_data[i][t] * (m.eta*m.mobility[i][j]/m.populations[i]) for j in m.mobility[i] if j in m.NODES)
        #return m.T_hat[i,t] == m.beta * (I_data[i]['I1'][t] + I_data[i]['I2'][t] + I_data[i]['I3'][t]) * S_data[i][t] / m.populations[i]

    model.infection_process = pe.Constraint(model.NODES, model.TIMES, rule=_infection_process)

    # least squares objective function
    def _lse(m, i):
        return sum( (m.T_hat[i,t] - m.T_data[i][t])**2 for t in m.TIMES)
    model.lse = pe.Expression(model.NODES, rule=_lse)

    def _total_lse(m):
        return sum( m.lse[i] for i in m.NODES )
    model.total_lse = pe.Objective(rule=_total_lse)

    ## likelihood objective function
    #def _like(m):
    #    #return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one if I1[t-1] + I2[t-1] + I3[t-1] > 0) + \
    #        #    sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one)
    #    return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one if I1[t] + I2[t] + I3[t] > 0) + \
    #        sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one)
    #model.like = pe.Objective(rule=_like, sense=pe.maximize)
    #model.like.deactivate()

    return model

