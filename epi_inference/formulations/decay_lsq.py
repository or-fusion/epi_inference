__all__ = ['run_decay_lsq']

import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import check_optimal_termination

def run_decay_lsq(*, recon, analysis_window, select_window=None, verbose=False):
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
    select_window: int
       An integer indicating that only this window will be used in the formulation.
    verbose : bool
       If true, then more output is printed to the console when the analysis is run
    """
    # create the Pyomo optimization formulation
    m = create_inference_formulation_decay(
        recon=recon,
        analysis_window=analysis_window,
        select_window=select_window,
        verbose=verbose
    )

    if m is None:
        return {'beta': None, 'status': 'failed', 'msg': 'Empty model.'}

    # call the solver
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m)

    # Check that the solve completed successfully
    if check_optimal_termination(status) == False:
        return {'beta': None,
                'status': 'failed',
                'msg': 'Unknown solver error.'
                }

    # check that the beta value was successfully solved
    if m.beta.stale == True:
        return {'est_beta': None,
                'status': 'stale',
                'msg': 'Transmission parameter beta not solved (stale).',
                }

    return {'est_beta': value(m.beta),
            'status': 'ok',
            'msg': 'Optimal solution found',
            }

def create_inference_formulation_decay(*, recon, analysis_window, select_window, verbose=False):
    """
    Creates a one-step-ahead inference model using a decay model with
    3 I compartments. The model is written in terms of absolute numbers
    of cases (not ln-transform)

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
    #model.DATES = dates
    #model.T_data = T
    #model.S_data = S
    #model.E_data = E
    #model.I1_data = I1
    #model.I2_data = I2
    #model.I3_data = I3
    #model.R_data = R

    # cache some parameters on the model
    #model.N = population
    #model.sigma = sigma
    #model.gamma_1 = gamma_1
    #model.gamma_2 = gamma_2
    #model.gamma_3 = gamma_3
    #model.report_delay = report_delay
    #model.rep_fac = reporting_factor

    # define the set of times
    model.timesteps = [i for i in range(len(T))]
    model = pe.ConcreteModel()
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
