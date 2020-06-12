import pandas as pd
from datetime import datetime
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination
from .decay_lsq import create_inference_formulation_decay


# ToDo: add datetime handling to pass in the dates associated with the data
def run_decay_blike(cm_rep_cases, population, sigma, gamma, report_delay, reporting_factor, analysis_window, Cdates, verbose=False):
    """
    This function solves a likelihood inference formulation
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
    assert sigma > 0
    assert gamma > 0
    assert report_delay > 0
    assert population > 0
    assert reporting_factor >= 1

    # initialize the likelihood model with results from the lsq model
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
    
    # solve the least-squares problem
    solver = SolverFactory('ipopt')
    status = solver.solve(m, tee=verbose)
    
    # Check that the initialization solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error in least-squares initialization.'}

    # check that the beta value was successfully solved
    if m.beta.stale == True:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Transmission parameter beta not solved (stale) in least-squares initialization.'}

    # deactivate the least-squares objective and activate the likelihood objective
    m.o_like.activate()
    m.o_lse.deactivate()

    # solve the likelihood formulation using initialization from the least-squares
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Unknown solver error.', 'population':population}

    # check that the beta value was successfully solved
    if m.beta.stale == True:
        return {'est_beta': None, 'status': 'failed', 'msg': 'Transmission parameter beta not solved (stale).', 'population':population}

    return {'est_beta': value(m.beta),'status': 'ok', 'msg': 'Optimal solution found', 'population':population}
