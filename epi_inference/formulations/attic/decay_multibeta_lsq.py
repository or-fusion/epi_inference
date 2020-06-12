import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import assert_optimal_termination
import math

from ..tseir_utils import compute_compartments_decay


def run_decay_multibeta_lsq(cm_rep_cases, population, sigma, gamma, deltaP, reporting_factor, verbose=False):
    print("WARNING: THIS MODEL IS NOT TESTED")
    assert sigma > 0
    assert gamma > 0
    assert deltaP > 0
    assert population > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_decay_multibeta(
            cumulative_reported_cases=cm_rep_cases,
            population=population, sigma=sigma,
            gamma_1=gamma*3, gamma_2=gamma*3, gamma_3=gamma*3,
            deltaP=deltaP, reporting_factor=reporting_factor,
            verbose=verbose,
            lse=True
    )

    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)
    assert_optimal_termination(status)

    results = {}
    for i in m.beta:
        results['est_beta_week'+str(i)] = value(m.beta[i])
    return results

def create_inference_formulation_decay_multibeta(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor, verbose=False, lse=True):
    """
    Creates a nonlinear one-step-ahead inference model using a decay model with 3 I compartments. The model is written in terms of absolute numbers of cases (not ln-transform)

    Parameters
    ----------
    cumulative_reported_cases : list of *new* cases reported in each time period
    population : the population of the node being considered
    sigma : float
       the rate constant for cases leaving the E compartment (1/incubation period)
    gamma : float
       the rate constant for leaving the I compartment. This will be multiplied
       by 3 when applied to each of the three I compartments
    deltaP : int
       the number of days between when someone is infected and when
       they will become a reported case (This should only shift the data
       and not impact the inference results.)
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported

    """
    cases, S, E, I1, I2, I3, R = compute_compartments_decay(cumulative_reported_cases=cumulative_reported_cases,
                                                            population=population, sigma=sigma,
                                                            gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                                            deltaP=deltaP, reporting_factor=reporting_factor)
    assert(len(cumulative_reported_cases) == len(cases))

    if verbose:
        print('corrected case data being used:')
        print(cases)
    
    model = pe.ConcreteModel()

    # cache some parameters on the model to make
    # reporting easier
    model.N = population # overall population
    model.sigma = sigma
    model.gamma_1 = gamma_1
    model.gamma_2 = gamma_2
    model.gamma_3 = gamma_3
    model.deltaP = deltaP
    model.rep_fac = reporting_factor

    # define the set of times
    model.timesteps = [i for i in range(len(cases))]
    model.TIMES = pe.Set(initialize=model.timesteps, ordered=True)
    model.TIMES_m_one = pe.Set(initialize=model.timesteps[1:], ordered=True)

    numcases = len(cumulative_reported_cases)
    model.TIMES_beta = pe.Set(initialize=[i for i in range((numcases-1+6)//7)], ordered=True)
    model.beta_offset = 7*((numcases-1+6)//7) - numcases + 1
    model.beta = pe.Var(model.TIMES_beta, initialize=1.3, bounds=(0,None)) # transmission parameter
    ##model.alpha = pe.Var(initialize=1.0)
    ##model.alpha.fix(1.0)

    counter = [0]*((numcases-1+6)//7)

    #print(len(cases))
    #print(cases)
    #print("numcases", numcases)
    #model.TIMES.pprint()
    #model.TIMES_beta.pprint()
    #print("offset", model.beta_offset)

    # define the case count variables
    model.Chat = pe.Var(model.TIMES_m_one, initialize=1.0)

    # infection process
    def _infection_process(m, t):
        if t == m.TIMES.last():
            return pe.Constraint.Skip
        counter[(t+model.beta_offset)//7] = counter[(t+model.beta_offset)//7] + 1
        return m.Chat[t+1] == model.beta[(t+model.beta_offset)//7] * (I1[t] + I2[t] + I3[t]) * S[t] / m.N
    model.infection_process = pe.Constraint(model.TIMES, rule=_infection_process)

    if lse:
        # least squares objective function
        def _lse(m):
            return sum( (m.Chat[t] - cases[t])**2 for t in m.TIMES_m_one)
        model.o_lse = pe.Objective(rule=_lse)

    else:
        # likelihood objective function

        #
        # Alternate likelihood function
        #
        #def _like(m):
            #return sum( cases[t]/m.N*pe.log(1-pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / m.N)) for t in m.TIMES_m_one if I1[t-1] + I2[t-1] + I3[t-1] > 0) + \
                #    sum( (S[t]-cases[t])/m.N*pe.log(pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / m.N)) for t in m.TIMES_m_one)

        def _like(m):
            return sum( cases[t]/m.N*pe.log(1-pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / m.N)) for t in m.TIMES_m_one if I1[t] + I2[t] + I3[t] > 0) + \
                sum( (S[t]-cases[t])/m.N*pe.log(pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / m.N)) for t in m.TIMES_m_one)
        model.o_like = pe.Objective(rule=_like, sense=pe.maximize)

    if verbose:
        print(counter)

    return model
