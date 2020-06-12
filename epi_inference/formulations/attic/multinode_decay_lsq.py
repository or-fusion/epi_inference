import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
import pandas as pd
from datetime import datetime
from . import reconstruction as recon

def run_multinode_decay_lsq(cm_rep_cases, populations, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
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

    # create the Pyomo optimization formulation
    m = create_inference_formulation_multinode_decay_lsq(
        Cdates=Cdates,
        cumulative_reported_cases=cm_rep_cases,
        populations=populations,
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

    # check that the beta value was successfully solved
    if m.beta.stale == True:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Transmission parameter beta not solved (stale).'}

    return {'est_beta': value(m.beta),'status': 'ok', 'msg': 'Optimal solution found'}

def create_inference_formulation_multinode_decay_lsq(Cdates, cumulative_reported_cases, populations, sigma, gamma_1, gamma_2, gamma_3, report_delay, reporting_factor, analysis_window, verbose=False):
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
    model.report_delay = report_delay
    model.reporting_factor = reporting_factor

    model.beta = pe.Var(initialize=1.3, bounds=(0,None)) # transmission parameter
    # for now, alpha is not used
    # model.alpha = pe.Var(initialize=1.0)
    # model.alpha.fix(1.0)

    #model.A = pe.Set(initialize=[v for v in cumulative_reported_cases.keys().to_list()])
    model.A = pe.Set(initialize=list(range(len(cumulative_reported_cases.keys()))))

    def block_rule(B, nodeid):
        # Cached data
        B.N = float(populations.iloc[nodeid]) # overall population
        cm_rep_cases_node = [v[0] for v in cumulative_reported_cases.iloc[:, [nodeid] ].values]

        rdates, rcases, dates, T, S, E, I1, I2, I3, R = \
            recon.reconstruct_states_deterministic_decay(
                Cdates=Cdates,
                cumulative_reported_cases=cm_rep_cases_node,
                population=B.N,
                sigma=sigma,
                gamma=gamma_1/3,
                reporting_factor=reporting_factor,
                report_delay=report_delay
            )

        if verbose:                                     # pragma: no cover
            print('corrected case data being used:')
            print(T)

        # define the set of times
        B.timesteps = [i for i in range(len(T))]
        B.TIMES = pe.Set(initialize=B.timesteps, ordered=True)
        #B.TIMES_m_one = pe.Set(initialize=B.timesteps[1:], ordered=True)

        # define the case count variables
        B.T_hat = pe.Var(B.TIMES, initialize=1.0)

        # infection process
        def _infection_process(b, t):
            return b.T_hat[t] == model.beta * (I1[t] + I2[t] + I3[t]) * S[t] / b.N
        B.infection_process = pe.Constraint(B.TIMES, rule=_infection_process)

        # least squares objective function
        def _lse(b):
            return sum( (b.T_hat[t] - T[t])**2 for t in b.TIMES)
        B.lse = pe.Expression(rule=_lse)

    model.b = pe.Block(model.A, rule=block_rule)

    def _total_lse(m):
        return sum( m.b[a].lse for a in m.A )
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

