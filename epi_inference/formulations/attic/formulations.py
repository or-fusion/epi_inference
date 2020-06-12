#
# THIS IS THE OLD FILE WITH FORMULATIONS
#
# These models have been moved to seperate files, except for the *delay*
# formulations, which are not actively used.
#
import pyomo.environ as pe
from pyomo.environ import SolverFactory, value
from pyomo.opt import assert_optimal_termination
import math

from ..tseir_utils import compute_compartments_time_delay, compute_compartments_decay


def run_multinode_decay_lsq(cm_rep_cases, population, sigma, gamma, deltaP, reporting_factor, verbose=False):
    assert sigma > 0
    assert gamma > 0
    assert deltaP > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_multinode_decay_lsq(
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
    #m.pprint()
    #m.display()

    return {'est_beta': value(m.beta)}

def run_decay_lsq(cm_rep_cases, population, sigma, gamma, deltaP, reporting_factor, verbose=False):
    assert sigma > 0
    assert gamma > 0
    assert deltaP > 0
    assert population > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_decay_lsq(
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

    return {'est_beta': value(m.beta)}

def run_decay_blike(cm_rep_cases, population, sigma, gamma, deltaP, reporting_factor, verbose=False):
    assert sigma > 0
    assert gamma > 0
    assert deltaP > 0
    assert population > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_decay_lsq(
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

    m.lse.deactivate()
    m.like.activate()
    solver = SolverFactory('ipopt')
    solver.options['tol']=1e-8
    status = solver.solve(m, tee=verbose)
    assert_optimal_termination(status)   

    return {'est_beta': value(m.beta)}

def run_delay_lsq(cm_rep_cases, population, deltaE, deltaI, deltaP, reporting_factor, verbose=False):
    assert deltaE > 0
    assert deltaI > 0
    assert deltaP > 0
    assert population > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_delay_lsq(
            cumulative_reported_cases=cm_rep_cases,
            population=population,
            deltaE=deltaE,
            deltaI=deltaI,
            deltaP=deltaP,
            reporting_factor=reporting_factor,
            verbose=verbose
    )
    
    solver = SolverFactory('ipopt')
    status = solver.solve(m, tee=verbose)
    assert_optimal_termination(status)
    return {'est_beta': value(m.beta)}

def run_delay_ln_lsq(cm_rep_cases, population, deltaE, deltaI, deltaP, reporting_factor, verbose=False):
    assert deltaE > 0
    assert deltaI > 0
    assert deltaP > 0
    assert population > 0
    assert reporting_factor >= 1

    m = create_inference_formulation_delay_ln_lsq(
            cumulative_reported_cases=cm_rep_cases,
            population=population,
            deltaE=deltaE,
            deltaI=deltaI,
            deltaP=deltaP,
            reporting_factor=reporting_factor,
            verbose=verbose
    )
    
    solver = SolverFactory('ipopt')
    status = solver.solve(m, tee=verbose)
    assert_optimal_termination(status)
    return {'est_ln_beta': value(m.ln_beta),
            'est_beta': math.exp(value(m.ln_beta))}

def create_inference_formulation_delay_ln_lsq(cumulative_reported_cases, population, deltaE, deltaI, deltaP, reporting_factor, verbose=False):
    """Creates a one-step-ahead inference model based on time delays in
    each compartment. The model is written in terms of the
    ln-transform of the cases (not absolute)

    Parameters
    ----------
    cumulative_reported_cases : list of *new* cases reported in each time period
    population : the population of the node being considered
    deltaE : int
       the number of days in the exposed compartment
    deltaI : int
       the number of days in the infectious compartment
    deltaP : int
       the number of days between when someone is infected and when
       they will become a reported case (This should only shift the data
       and not impact the inference results.)
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported

    """
    cases, S, E, I, R = compute_compartments_time_delay(cumulative_reported_cases, population, deltaE, deltaI, deltaP, reporting_factor)

    if verbose:
        print('corrected case data being used:')
        print(cases)
    

    model = pe.ConcreteModel()

    # cache some parameters on the model to make
    # reporting easier
    model.N = population # overall population
    model.deltaE = deltaE
    model.deltaI = deltaI
    model.deltaP = deltaP
    model.rep_fac = reporting_factor

    # define the set of times
    model.timesteps = [i for i in range(len(cases))]
    model.TIMES = pe.Set(initialize=model.timesteps, ordered=True)
    model.TIMES_m_one = pe.Set(initialize=model.timesteps[1:], ordered=True)

    # define the parameters
    model.ln_beta = pe.Var(initialize=math.log(1.3)) # transmission parameter
    model.alpha = pe.Var(initialize=1.0)
    model.alpha.fix(1.0)

    # define the case count variables
    model.ln_Chat = pe.Var(model.TIMES_m_one, initialize=1.0)

    # log of the infection process
    def _infection_process(m, t):
        if t == m.TIMES.last():
            return pe.Constraint.Skip
        if I[t]==0:
            # if there are no infectives today, I am getting no information about beta
            # with this equation - let ln_Chat be free
            if verbose:
                print(' *** No infectives at time ', t, 'skipping equation for ln_Chat[', t+1, ']')
            return pe.Constraint.Skip
        return m.ln_Chat[t+1] == m.ln_beta + m.alpha*math.log(I[t]) + math.log(S[t]) - math.log(model.N)
    model.infection_process = pe.Constraint(model.TIMES, rule=_infection_process)

    # least squares objective function
    def _lse(m):
        """
        expr = 0
        for t in model.TIMES_m_one:
            if cases[t-1] > 0:
                # we have computed an ln_Chat at t
                if cases[t] == 0:
                    # we computed a ln_Chat at t, but there are no cases in
                    # the next timestep. Therefore, match a low number
                    expr += (model.ln_Chat[t] - math.log(1e-8))**2
                else:
                    # we computed a ln_Chat at t, and there are cases in
                    # the next timestep.                    
                    expr += (model.ln_Chat[t] - math.log(cases[t]))**2
        return expr
        """
 
        expr = 0
        for t in model.TIMES_m_one:
            if cases[t] > 0:
                expr += (model.ln_Chat[t] - math.log(cases[t]))**2
            elif I[t-1] > 0: # There were cases to compute ln_Chat[t], but no cases now
                # we had cases before, but we do not now - need to think about the
                # value to include in the log below
                #print('TIMESTEP', t)
                expr += (model.ln_Chat[t] - math.log(0.75))**2
        return expr
 
        # TODO: Look above - instead of skipping, we may need to "match beta"
#        return sum( (model.ln_Chat[t] - math.log(cases[t]))**2 for t in model.TIMES_m_one if cases[t] > 0)
#        return sum( (model.ln_Chat[t] - math.log(cases[t]))**2 for t in model.TIMES_m_one)
    model.lse = pe.Objective(rule=_lse)
    model.lse.pprint()
    
    return model


def create_inference_formulation_delay_lsq(cumulative_reported_cases, population, deltaE, deltaI, deltaP, reporting_factor, verbose=False):
    """Creates a one-step-ahead inference model based on time delays in
    each compartment. The model is written in terms of the
    absolute number of cases (not ln-transform)

    Parameters
    ----------
    cumulative_reported_cases : list of *new* cases reported in each time period
    population : the population of the node being considered
    deltaE : int
       the number of days in the exposed compartment
    deltaI : int
       the number of days in the infectious compartment
    deltaP : int
       the number of days between when someone is infected and when
       they will become a reported case (This should only shift the data
       and not impact the inference results.)
    reporting_factor : float
       The reporting factor (>1). If set to 5 this means 1 in 5 cases is reported
    """
    cases, S, E, I, R = compute_compartments_time_delay(cumulative_reported_cases, population, deltaE, deltaI, deltaP, reporting_factor)

    if verbose:
        print('corrected case data being used:')
        print(cases)
    
    model = pe.ConcreteModel()

    # cache some parameters on the model to make
    # reporting easier
    model.N = population # overall population
    model.deltaE = deltaE
    model.deltaI = deltaI
    model.deltaP = deltaP
    model.rep_fac = reporting_factor

    # define the set of times
    model.timesteps = [i for i in range(len(cases))]
    model.TIMES = pe.Set(initialize=model.timesteps, ordered=True)
    model.TIMES_m_one = pe.Set(initialize=model.timesteps[1:], ordered=True)

    # define the parameters
    model.beta = pe.Var(initialize=1.3) # transmission parameter
    model.alpha = pe.Var(initialize=1.0)
    model.alpha.fix(1.0)

    # define the case count variables
    model.Chat = pe.Var(model.TIMES_m_one, initialize=1.0)

    # log of the infection process
    def _infection_process(m, t):
        if t == m.TIMES.last():
            return pe.Constraint.Skip
        return m.Chat[t+1] == m.beta * I[t]*S[t]/model.N
    model.infection_process = pe.Constraint(model.TIMES, rule=_infection_process)

    # least squares objective function
    def _lse(m):
        return sum( (model.Chat[t] - cases[t])**2 for t in model.TIMES_m_one)
#        return sum( 1.0/max(cases[t]**2,1.0)*(model.Chat[t] - cases[t])**2 for t in model.TIMES_m_one)
    model.lse = pe.Objective(rule=_lse)

    return model


def create_inference_formulation_multinode_decay_lsq(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor, verbose=False):
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

    model.beta = pe.Var(initialize=1.3, bounds=(0,None)) # transmission parameter
    ## IS THIS AN ERROR?
    ##model.alpha = pe.Var(initialize=1.0)
    ##model.alpha.fix(1.0)

    model.A = pe.Set(initialize=[v for v in cumulative_reported_cases.keys().to_list()])

    def block_rule(B, nodeid):
        # Cached data
        B.N = population[int(nodeid)] # overall population

        cases, S, E, I1, I2, I3, R = compute_compartments_decay(cumulative_reported_cases=cumulative_reported_cases[nodeid],
                                                            population=B.N,
                                                            sigma=sigma,
                                                            gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                                            deltaP=deltaP, reporting_factor=reporting_factor)
        if verbose:
            print('corrected case data being used:')
            print(cases)

        # define the set of times
        B.timesteps = [i for i in range(len(cases))]
        B.TIMES = pe.Set(initialize=B.timesteps, ordered=True)
        B.TIMES_m_one = pe.Set(initialize=B.timesteps[1:], ordered=True)

        # define the case count variables
        B.Chat = pe.Var(B.TIMES_m_one, initialize=1.0)

        # infection process
        def _infection_process(b, t):
            if t == b.TIMES.last():
                return pe.Constraint.Skip
            return b.Chat[t+1] == model.beta * (I1[t] + I2[t] + I3[t]) * S[t] / b.N
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

    return model

def create_inference_formulation_decay_lsq(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor, verbose=False):
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

    # define the parameters
    model.beta = pe.Var(initialize=1.3, bounds=(0,None)) # transmission parameter
    model.alpha = pe.Var(initialize=1.0)
    model.alpha.fix(1.0)

    # define the case count variables
    model.Chat = pe.Var(model.TIMES_m_one, initialize=1.0)

    # infection process
    def _infection_process(m, t):
        if t == m.TIMES.last():
            return pe.Constraint.Skip
        return m.Chat[t+1] == m.beta * (I1[t] + I2[t] + I3[t]) * S[t] / model.N
    model.infection_process = pe.Constraint(model.TIMES, rule=_infection_process)

    # least squares objective function
    def _lse(m):
        return sum( (model.Chat[t] - cases[t])**2 for t in model.TIMES_m_one)
    model.lse = pe.Objective(rule=_lse)

    # likelihood objective function
    def _like(m):
        #return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one if I1[t-1] + I2[t-1] + I3[t-1] > 0) + \
            #    sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one)
        return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one if I1[t] + I2[t] + I3[t] > 0) + \
            sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one)
    model.like = pe.Objective(rule=_like, sense=pe.maximize)
    model.like.deactivate()

    return model

def create_inference_formulation_decay_blike(cumulative_reported_cases, population, sigma, gamma_1, gamma_2, gamma_3, deltaP, reporting_factor, verbose=False):
    """
    Creates a nonlinear one-step-ahead inference model using a decay model with 3 I compartments. The model is written in terms of absolute numbers of cases (not ln-transform).
    This formulation uses a log-likelihood for the binomial distribution to estimate beta.

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

    # define the parameters
    model.beta = pe.Var(initialize=1.3, bounds=(0,None)) # transmission parameter
    model.alpha = pe.Var(initialize=1.0)
    model.alpha.fix(1.0)

    def _like(m):
        #return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one if I1[t-1] + I2[t-1] + I3[t-1] > 0) + \
        #    sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t-1] + I2[t-1] + I3[t-1]) / model.N)) for t in model.TIMES_m_one)
        return sum( cases[t]/model.N*pe.log(1-pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one if I1[t] + I2[t] + I3[t] > 0) + \
            sum( (S[t]-cases[t])/model.N*pe.log(pe.exp(-m.beta * (I1[t] + I2[t] + I3[t]) / model.N)) for t in model.TIMES_m_one)
    model.like = pe.Objective(rule=_like, sense=pe.maximize)
    
    return model
