import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import assert_optimal_termination
import math

from ..tseir_utils import compute_compartments_decay


def run_multinode_decay_multibeta_lsq(cm_rep_cases, population, sigma, gamma, deltaP, reporting_factor, verbose=False):
    assert sigma > 0
    assert gamma > 0
    assert deltaP > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_multinode_decay_multibeta_lsq(
            cumulative_reported_cases=cm_rep_cases,
            population=population, sigma=sigma,
            gamma_1=gamma*3, gamma_2=gamma*3, gamma_3=gamma*3,
            deltaP=deltaP, reporting_factor=reporting_factor,
            verbose=verbose
    )

    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)
    assert_optimal_termination(status)

    results = {}
    for i in m.beta:
        results['est_beta_week'+str(i)] = value(m.beta[i])
    return results

def create_inference_formulation_multinode_decay_multibeta_lsq(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor, verbose=False):
    """
    Creates a nonlinear one-step-ahead inference model using a decay
    model with 3 I compartments. The model is written in terms of absolute
    numbers of cases (not ln-transform).  The model combines estimates across
    multiple time series, one for each node.

    Parameters
    ----------

    cumulative_reported_cases : a dataframe of *new* cases reported in
       each time period; each column in the dataframe is a separate time
       series
    population : a dataframe with a single column that represents the
       population for different columns in cumulative_reported_cases
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
    model = pe.ConcreteModel()

    # Cached data
    model.sigma = sigma
    model.gamma_1 = gamma_1
    model.gamma_2 = gamma_2
    model.gamma_3 = gamma_3
    model.deltaP = deltaP
    model.reporting_factor = reporting_factor

    #
    # NOTE: Here we assume that all time series have the same length
    #
    numcases = cumulative_reported_cases.shape[0]
    model.TIMES_beta = pe.Set(initialize=[i for i in range((numcases-1+6)//7)], ordered=True)
    model.beta_offset = 7*((numcases-1+6)//7) - numcases + 1
    model.beta = pe.Var(model.TIMES_beta, initialize=1.3, bounds=(0,None)) # transmission parameter
    ## IS THIS AN ERROR?
    ##model.alpha = pe.Var(initialize=1.0)
    ##model.alpha.fix(1.0)

    counter = [0]*((numcases-1+6)//7)
    #print("len(cases)", len(cumulative_reported_cases['12121']))
    #print("len(cases)", len(cumulative_reported_cases['12121'].to_list()))
    #print("cases", cumulative_reported_cases['12121'].to_list())
    #print("HERE", numcases, 7*((numcases+6)//7), model.beta_offset)
    #model.TIMES_beta.pprint()
    #print(counter)

    #model.A = pe.Set(initialize=[v for v in cumulative_reported_cases.keys().to_list()])
    model.A = pe.Set(initialize=list(range(len(cumulative_reported_cases.keys()))))

    def block_rule(B, nodeid):
        # Cached data
        B.N = population.iloc[nodeid] # overall population
        reported_cases = [v[0] for v in cumulative_reported_cases.iloc[:, [nodeid] ].values]

        cases, S, E, I1, I2, I3, R = compute_compartments_decay(cumulative_reported_cases=reported_cases,
                                                            population=B.N,
                                                            sigma=sigma,
                                                            gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                                            deltaP=deltaP, reporting_factor=reporting_factor)
        assert(len(cases) == len(reported_cases))

        if verbose:
            print('corrected case data being used:')
            print(cases)

        # define the set of times
        B.timesteps = [i for i in range(len(cases))]
        B.TIMES = pe.Set(initialize=B.timesteps, ordered=True)
        B.TIMES_m_one = pe.Set(initialize=B.timesteps[1:], ordered=True)

        #print("id", nodeid)
        #print("LEN(CASES)", len(cases))
        #print("cases", cases)
        #print("TIMES")
        #B.TIMES.pprint()

        # define the case count variables
        B.Chat = pe.Var(B.TIMES_m_one, initialize=1.0)

        # infection process
        def _infection_process(b, t):
            if t == b.TIMES.last():
                return pe.Constraint.Skip
            counter[(t+model.beta_offset)//7] = counter[(t+model.beta_offset)//7] + 1
            return b.Chat[t+1] == model.beta[(t+model.beta_offset)//7] * (I1[t] + I2[t] + I3[t]) * S[t] / b.N
        B.infection_process = pe.Constraint(B.TIMES, rule=_infection_process)

        # least squares objective function
        def _lse(b):
            return sum( (b.Chat[t] - cases[t])**2 for t in b.TIMES_m_one)
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

    if verbose:
        print("counter",counter)

    return model

