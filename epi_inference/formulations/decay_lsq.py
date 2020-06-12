import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
import pandas as pd
from datetime import datetime
from epi_inference.reconstruction import common as rcommon
from epi_inference.reconstruction import deterministic as recond

# ToDo: add datetime handling to pass in the dates associated with the data
def run_decay_lsq(cm_rep_cases, population, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
    """
    This function solves the least-squares inference inference formulation
    using the decay-based reconstruction function.

    Parameters
    ----------
    cm_rep_cases : list of *new* cases reported in each time period
       Note that this list is 1 entry longer than the transmissions, and 
       it must start at zero (based on the assumptions in the reconstruction).
    population : the population of the node being considered
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
    # check validity of some of the inputs
    assert sigma > 0
    assert gamma > 0
    assert report_delay > 0
    assert population > 0
    assert reporting_factor >= 1

    # create the Pyomo optimization formulation
    m = create_inference_formulation_decay(
        Cdates=Cdates,
        cumulative_reported_cases=cm_rep_cases,
        population=population, sigma=sigma,
        gamma_1=gamma*3, gamma_2=gamma*3, gamma_3=gamma*3,
        report_delay=report_delay, reporting_factor=reporting_factor,
        analysis_window=analysis_window,
        verbose=verbose,
        lse=True
    )

    reconstruction = {
        'date': m.DATES,
        'S': m.S_data,
        'T': m.T_data,
        'E': m.E_data,
        'I1': m.I1_data,
        'I2': m.I2_data,
        'I3': m.I3_data,
        'R': m.R_data
    }
    
    # call the solver
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error.', 'population':population, 'total_cases':cm_rep_cases[-1], 'reconstruction': reconstruction}

    # check that the beta value was successfully solved
    if m.beta.stale == True:
        return {'est_beta': None,
                'status': 'stale',
                'msg': 'Transmission parameter beta not solved (stale).',
                'population':population,
                'total_cases':cm_rep_cases[-1],
                'reconstruction': reconstruction
        }

    return {'est_beta': value(m.beta),
            'status': 'ok',
            'msg': 'Optimal solution found',
            'population':population,
            'total_cases':cm_rep_cases[-1],
            'reconstruction': reconstruction
            }

def create_inference_formulation_decay(Cdates, cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, report_delay, reporting_factor, analysis_window, verbose=False, lse=True):
    """
    Creates a one-step-ahead inference model using a decay model with 3 I compartments. The model is written in terms of absolute numbers of cases (not ln-transform)

    Parameters
    ----------
    Cdates: list of datetime objects
       The list of datetime objects that correspond to the dates for the
       cumulative_reported_cases
    cumulative_reported_cases : list of *new* cases reported in each time period
    population : the population of the node being considered
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
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported
    analysis_window : dict or None
       This is a dictionary indicating the window of time that should be used 
       in the objective function. If None, then the full set of data will be used.
       The key "days" indicates the number of days from the end of the data that 
       should be used in the objective function.
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    lse : bool
       If true, the activated objective corresponds to the least-squares, otherwise the 
       likelihood objective will be activated.
    """
    rdates, rcases = rcommon.reported_cases_from_cumulative(Cdates, cumulative_reported_cases)
    dates, T, S, E, I1, I2, I3, R = recond.reconstruct_states_deterministic_decay(
        Cdates=Cdates,
        cumulative_reported_cases=cumulative_reported_cases,
        population=population,
        sigma=sigma,
        gamma=gamma_1/3,
        reporting_factor=reporting_factor,
        report_delay=report_delay
        )

    """
    T, S, E, I1, I2, I3, R = compute_compartments_decay(cumulative_reported_cases=cumulative_reported_cases,
                                                            population=population, sigma=sigma,
                                                            gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                                            report_delay=report_delay, reporting_factor=reporting_factor)
    """

    if verbose:                                     # pragma: no cover
        print('corrected case data being used:')
        print(T)
    
    model = pe.ConcreteModel()
    model.DATES = dates
    model.T_data = T
    model.S_data = S
    model.E_data = E
    model.I1_data = I1
    model.I2_data = I2
    model.I3_data = I3
    model.R_data = R

    # cache some parameters on the model
    model.N = population
    model.sigma = sigma
    model.gamma_1 = gamma_1
    model.gamma_2 = gamma_2
    model.gamma_3 = gamma_3
    model.report_delay = report_delay
    model.rep_fac = reporting_factor

    # define the set of times
    model.timesteps = [i for i in range(len(T))]
    model.TIMES = pe.Set(initialize=model.timesteps, ordered=True)

    start=0
    if analysis_window.get('days',None) is not None:
        start = max(start, len(model.timesteps)-analysis_window.get('days',None)-1)
    model.TIMES_m_obj = pe.Set(initialize=model.timesteps[start:], ordered=True)

    # define the parameters
    model.beta = pe.Var(initialize=1.3, bounds=(0,None)) # transmission parameter
    model.alpha = pe.Var(initialize=1.0)
    model.alpha.fix(1.0)

    # define the case count variables
    model.T_hat = pe.Var(model.TIMES, initialize=1.0)

    # infection process
    def _infection_process(m, t):
        return m.T_hat[t] == m.beta * (I1[t] + I2[t] + I3[t]) * S[t] / m.N
    model.infection_process = pe.Constraint(model.TIMES_m_obj, rule=_infection_process)

    # least squares objective function
    def _lse(m):
        return sum( (m.T_hat[t] - T[t])**2 for t in m.TIMES_m_obj)
    model.o_lse = pe.Objective(rule=_lse)

    def _like(m):
        return sum( T[t]/m.N*pe.log(1-pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / m.N)) for t in m.TIMES_m_obj if I1[t] + I2[t] + I3[t] > 0) + \
            sum( (S[t]-T[t])/m.N*pe.log(pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / m.N)) for t in m.TIMES_m_obj)
    model.o_like = pe.Objective(rule=_like, sense=pe.maximize)

    if lse:
        model.o_like.deactivate()
    else:
        model.o_lse.deactivate()
    return model
