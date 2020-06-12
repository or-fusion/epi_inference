import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
import pandas as pd
from datetime import datetime
from . import reconstruction as recon

def run_multibeta_multiwin_decay_lsq(cm_rep_cases, populations, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.  

    Each county has a separate beta value, but a common set of omega
    multipliers.

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
    m = create_inference_formulation_multibeta_multiwin_decay_lsq(
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

    #m.pprint()
    #m.display()

    est_beta = {}

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': est_beta, 'status': 'failed', 'msg': 'Unknown solver error.'}

    cm_rep_cases_node = [v[0] for v in cm_rep_cases.iloc[:, [0] ].values]

    rdates, rcases, dates, T, S, E, I1, I2, I3, R = \
        recon.reconstruct_states_deterministic_decay(
            Cdates=Cdates,
            cumulative_reported_cases=cm_rep_cases_node,
            population=populations.iloc[0],
            sigma=sigma,
            gamma=gamma/3,
            reporting_factor=reporting_factor,
            report_delay=report_delay
        )

    ##print(dates)
    ##print(Cdates)
    ##print("X", len(dates), len(Cdates))
    #
    # TODO - confirm that we should use dates[j] for the date index
    #
    for i in m.A:
        for j in m.b[i].beta:
            if not m.b[i].beta[j].stale:
                fips = populations.index[i]
                if not fips in est_beta:
                    est_beta[fips] = {}
                #print("HERE",j,dates[j])
                est_beta[fips][dates[j]] = value(m.b[i].beta[j])

    return {'est_beta': est_beta, 'status': 'ok', 'msg': 'Optimal solution found'}

def create_inference_formulation_multibeta_multiwin_decay_lsq(Cdates, cumulative_reported_cases, populations, sigma, gamma_1, gamma_2, gamma_3, report_delay, reporting_factor, analysis_window, verbose=False):
    """
    Creates a one-step-ahead inference model using a decay
    model with 3 I compartments. The model is written in terms of absolute
    numbers of cases (not ln-transform).  

    The model combines computes a different beta value for each time series, using a common set of 
    omega factors to describe the changes in beta over time for all counties.

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
    window = int(analysis_window.get('days',7))
    assert(window >= 1)

    model = pe.ConcreteModel()

    # Cached data
    model.sigma = sigma
    model.gamma_1 = gamma_1
    model.gamma_2 = gamma_2
    model.gamma_3 = gamma_3
    model.report_delay = report_delay
    model.reporting_factor = reporting_factor

    #
    # We assume that all counties have the same number of time steps
    #
    # TODO: This should be checked in inference.py
    #
    model.timesteps = [i for i in range(len(Cdates)-1)]
    model.TIMES = pe.Set(initialize=model.timesteps, ordered=True)

    #model.omega_timesteps = [i+window-1 for i in range(len(Cdates)-window)]
    #model.omega = pe.Var(model.omega_timesteps, initialize=1.3, bounds=(0,None)) # transmission parameter
    #model.omega_sum = pe.Constraint(expr=sum(model.omega[i] for i in model.omega) == len(model.omega))

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

        # Q: Is there a better way to do this?
        TIMES = model.TIMES

        # transmission parameter
        B.beta = pe.Var(TIMES, initialize=0.1, bounds=(0,None))
        B.lse_beta = pe.Expression(TIMES, initialize=0)

        first = None
        for i in model.TIMES:
            rhscoef =  (I1[i] + I2[i] + I3[i]) * S[i] / B.N 
            if rhscoef > 0:
                first = i
                break

        terms = []
        if not first is None:
            i = first + window-1
            while i < len(model.TIMES):
                #print(nodeid,first,i)
                for offset in range(window):
                    tau = i-offset
                    rhscoef =  (I1[tau] + I2[tau] + I3[tau]) * S[tau] / B.N 
                    terms.append((B.beta[i] * rhscoef - T[tau])**2)
                    #print(tau, (B.beta[i] * rhscoef - T[tau])**2)
                i = i + 1

        B.lse = pe.Expression(expr=sum(terms))

    model.b = pe.Block(model.A, rule=block_rule)

    def _total_lse(m):
        return sum( m.b[a].lse for a in m.A )
    model.total_lse = pe.Objective(rule=_total_lse)

    return model

