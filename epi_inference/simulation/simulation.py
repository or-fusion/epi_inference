"""
This module contains some simple simulation models for testing the 
inference methods
"""
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
from scipy.integrate import odeint
from pyutilib.misc.misc import Bunch

# the '*' forces the caller to use named keyword arguments
def simulate_continuous_seiiir_deterministic(*, y0, tf, beta, sigma, gamma, rho, N, report_delay, tx):
    """
    This method simulates a deterministic seiiir model based on the given parameters.
    The model is continuous time and uses a scipy integrator
    
    Parameters
    ----------
    y0 : dict
       Initial values for the compartments. Should include keys for 'S', 'E', 'I1' ,'I2', 'I3', 'R'
    tf : int
       The number of days to simulate
    beta : float
       Transmission rate in contacts per day
    sigma : float
       This is the rate of transfer out of the exposed compartment (e.g., 1/(incubation period))
    gamma : float
       This model includes 3 infectious compartments. Therefore this is an approximation of the 
       overall rate of transfer out of the I compartment. Here, we create three compartments,
       each with a rate of 3*gamma (where gamma is 1/(infectious period))
    rho : float
       This is the reporting factor. If rho=8, then 1 in 8 cases are reported
    N : number
       This is the overall population
    report_delay : int
       This is the number of days between infection and reporting
    tx : list of numbers or None
       The number of transmissions from external sources for each timestep. If set to None, then
       this is assumed to be zero
    """
    assert tf > 1 and type(tf) is int
    assert N > 1
    assert y0['S'] <= N
    assert y0['E'] <= N
    assert y0['I1'] <= N
    assert y0['I2'] <= N
    assert y0['I3'] <= N
    assert y0['R'] <= N
    assert y0['S'] >= 0
    assert y0['E'] >= 0
    assert y0['I1'] >= 0
    assert y0['I2'] >= 0
    assert y0['I3'] >= 0
    assert y0['R'] >= 0
    assert beta >= 0
    assert sigma >=0
    assert gamma >= 0
    assert rho >= 1
    assert report_delay >= 1 and type(report_delay) is int
    assert tf > report_delay+1
    assert tx is None # not currently supported

    def model(y, t, beta, sigma, gamma, N):
        S=y[0]
        E=y[1]
        I1=y[2]
        I2=y[3]
        I3=y[4]
        R=y[5]
        C=y[6]

        rhs = np.zeros(7)
        rhs[0] = -beta*(I1+I2+I3)*S/N
        rhs[1] = beta*(I1+I2+I3)*S/N - sigma*E
        rhs[2] = sigma*E - 3*gamma*I1
        rhs[3] = 3*gamma*I1 - 3*gamma*I2
        rhs[4] = 3*gamma*I2 - 3*gamma*I3
        rhs[5] = 3*gamma*I3
        rhs[6] = beta*(I1+I2+I3)*S/N
        return rhs

    y = np.zeros(7)
    y[0] = y0['S']
    y[1] = y0['E']
    y[2] = y0['I1']
    y[3] = y0['I2']
    y[4] = y0['I3']
    y[5] = y0['R']
    y[6] = 0
    
    times = list(range(tf))

    y = odeint(lambda y,t : model(y, t, beta, sigma, gamma, N), y, times)
    S = y[:-1,0].tolist()
    E = y[:-1,1].tolist()
    I1 = y[:-1,2].tolist()
    I2 = y[:-1,3].tolist()
    I3 = y[:-1,4].tolist()
    R = y[:-1,5].tolist()
    C = y[:,6].tolist()
    T = [C[i+1]-C[i] for i in range(len(C)-1)]
    
    Cdates = pd.date_range(end=datetime(year=2020, month=5, day=15), periods=len(C)).to_pydatetime().tolist()
    dates = list()
    delta = timedelta(days=report_delay-1) # the -1 is because we add a 0 to the beginning of cm_rep_cases (to start at 0)
    for i in range(len(Cdates)-1):
        dates.append(Cdates[i]-delta)

    cumulative_reported_cases = Bunch(dates=Cdates, values=C)
    SEIIIR = Bunch(dates=dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=T)
    return Bunch(cumulative_reported_cases=cumulative_reported_cases, SEIIIR=SEIIIR)
    

def simulate_discrete_seiiir_deterministic(y0, tf, beta, sigma, gamma, rho, N, report_delay, tx):
    """
    This method simulates a deterministic seiiir model based on the given parameters.
    The model is discrete time and discretized by days
    
    Parameters
    ----------
    y0 : dict
       Initial values for the compartments. Should include keys for 'S', 'E', 'I1' ,'I2', 'I3', 'R'
    tf : int
       The number of days to simulate
    beta : float
       Transmission rate in contacts per day
    sigma : float
       This is the rate of transfer out of the exposed compartment (e.g., 1/(incubation period))
    gamma : float
       This model includes 3 infectious compartments. Therefore this is an approximation of the 
       overall rate of transfer out of the I compartment. Here, we create three compartments,
       each with a rate of 3*gamma (where gamma is 1/(infectious period))
    rho : float
       This is the reporting factor. If rho=8, then 1 in 8 cases are reported
    N : number
       This is the overall population
    report_delay : int
       This is the number of days between infection and reporting
    tx : list of numbers or None
       The number of transmissions from external sources for each timestep. If set to None, then
       this is assumed to be zero
    """
    assert tf > 1 and type(tf) is int
    assert N > 1
    assert y0['S'] <= N
    assert y0['E'] <= N
    assert y0['I1'] <= N
    assert y0['I2'] <= N
    assert y0['I3'] <= N
    assert y0['R'] <= N
    assert y0['S'] >= 0
    assert y0['E'] >= 0
    assert y0['I1'] >= 0
    assert y0['I2'] >= 0
    assert y0['I3'] >= 0
    assert y0['R'] >= 0
    assert beta >= 0
    assert sigma >=0
    assert gamma >= 0
    assert rho >= 1
    assert report_delay >= 1 and type(report_delay) is int
    assert tf > report_delay+1
    if tx is None:
        tx = [0]*tf

    S = [None]*tf
    E = [None]*tf
    T = [None]*tf # new transmissions
    I1 = [None]*tf
    I2 = [None]*tf
    I3 = [None]*tf
    R = [None]*tf

    S[0] = y0['S']
    E[0] = y0['E']
    I1[0] = y0['I1']
    I2[0] = y0['I2']
    I3[0] = y0['I3']
    R[0] = y0['R']

    cm_rep_cases = [None]*(tf+1)
    cm_rep_cases[0] = 0

    T[0] = beta * (I1[0] + I2[0] + I3[0]) * S[0] / N + tx[0]
    for t in range(0,tf-1):
        #print(T[t], S[t], E[t], I1[t], I2[t],I3[t],cm_rep_cases[t])
        assert T[t] >= 0
        S[t+1] = S[t] - T[t]
        S[t+1] = max(0,S[t+1])

        E[t+1] = E[t] + T[t] - sigma*E[t]
        E[t+1] = max(0, E[t+1])

        I1[t+1] = I1[t] + sigma*E[t] - gamma*3*I1[t]
        I1[t+1] = max(0, I1[t+1])

        I2[t+1] = I2[t] + gamma*3*I1[t] - gamma*3*I2[t]
        I2[t+1] = max(0, I2[t+1])

        I3[t+1] = I3[t] + gamma*3*I2[t] - gamma*3*I3[t]
        I3[t+1] = max(0, I3[t+1])

        R[t+1] = R[t] + gamma*3*I3[t]
        assert R[t+1] >= 0

        cm_rep_cases[t+1] = cm_rep_cases[t] + 1/rho*T[t]
        T[t+1] = beta * (I1[t+1] + I2[t+1] + I3[t+1]) * S[t+1] / N + tx[t+1]

    cm_rep_cases[tf] = cm_rep_cases[tf-1] + 1/rho*T[tf-1]
    cm_rep_cases_dates = pd.date_range(end=datetime(year=2020, month=5, day=15), periods=len(cm_rep_cases)).to_pydatetime().tolist()
    dates = list()
    delta = timedelta(days=report_delay-1) # the -1 is because we add a 0 to the beginning of cm_rep_cases (to start at 0)
    for i in range(len(cm_rep_cases_dates)-1):
        dates.append(cm_rep_cases_dates[i]-delta)

    cumulative_reported_cases = Bunch(dates=cm_rep_cases_dates, values=cm_rep_cases)
    SEIIIR = Bunch(dates=dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=T)
    return Bunch(cumulative_reported_cases=cumulative_reported_cases, SEIIIR=SEIIIR)


def simulate_discrete_seiiir_stochastic(y0, tf, beta, sigma, gamma, rho, N, report_delay, tx):
    assert tf > 1 and type(tf) is int
    assert N > 1
    assert y0['S'] <= N
    assert y0['E'] <= N
    assert y0['I1'] <= N
    assert y0['I2'] <= N
    assert y0['I3'] <= N
    assert y0['R'] <= N
    assert y0['S'] >= 0
    assert y0['E'] >= 0
    assert y0['I1'] >= 0
    assert y0['I2'] >= 0
    assert y0['I3'] >= 0
    assert y0['R'] >= 0
    assert beta >= 0
    assert sigma >=0
    assert gamma >= 0
    assert rho >= 1
    assert report_delay >= 1 and type(report_delay) is int
    assert tf > report_delay+1
    if tx is None:
        tx = [0]*tf


    # expand the parameters if necessary
    assert type(beta) is float or (type(beta) is np.ndarray and beta.dtype == np.float64)
    assert type(sigma) is float or (type(sigma) is np.ndarray and sigma.dtype == np.float64)
    assert type(gamma) is float or (type(gamma) is np.ndarray and gamma.dtype == np.float64)
    #assert type(report_delay) is int or (type(report_delay) is np.ndarray and report_delay.dtaype == np.int64)
    #assert type(rho) is float or (type(rho) is np.ndarray and rho.dtype == np.int64)

    if type(beta) is float:
        beta = beta*np.ones(tf)
    if type(sigma) is float:
        sigma = sigma*np.ones(tf)
    if type(gamma) is float:
        gamma = gamma*np.ones(tf)
    #if type(report_delay) is int:
    #    report_delay = report_delay*np.ones(tf).astype(np.int64)
    if type(rho) is float:
        rho = rho*np.ones(tg)

    S = [None]*tf
    E = [None]*tf
    T = [None]*tf # new transmissions
    I1 = [None]*tf
    I2 = [None]*tf
    I3 = [None]*tf
    R = [None]*tf

    S[0] = y0['S']
    E[0] = y0['E']
    I1[0] = y0['I1']
    I2[0] = y0['I2']
    I3[0] = y0['I3']
    R[0] = y0['R']

    p_EI = 1-np.exp(-sigma)
    p_II = 1-np.exp(-3*gamma)
    for t in range(0,tf-1):
        p_expose = 1-np.exp(-beta[t]*(I1[t]+I2[t]+I3[t])/N)
        T[t] = np.random.binomial(S[t], p_expose) + tx[t]
        Etout = np.random.binomial(E[t], p_EI[t])
        I1tout = np.random.binomial(I1[t], p_II[t])
        I2tout = np.random.binomial(I2[t], p_II[t])
        I3tout = np.random.binomial(I3[t], p_II[t])
        S[t+1] = S[t] - T[t]
        E[t+1] = E[t] + T[t] - Etout
        I1[t+1] = I1[t] + Etout - I1tout
        I2[t+1] = I2[t] + I1tout - I2tout
        I3[t+1] = I3[t] + I2tout - I3tout
        R[t+1] = R[t] + I3tout

    # get the last entry for transmissions
    T[tf-1] = np.random.binomial(S[tf-1], p_expose) + tx[tf-1]

    cm_rep_cases = [0]*(tf+1)
    for tp in range(1,tf+1):
        cm_rep_cases[tp] = cm_rep_cases[tp-1] + 1/rho*T[tp-1]

    # create and shift the dates to account for reporting delay
    cm_rep_cases_dates = pd.date_range(end=datetime(year=2020, month=5, day=15), periods=len(cm_rep_cases)).to_pydatetime().tolist()
    dates = list()
    delta = timedelta(days=report_delay-1) # the -1 is because we add a 0 to the beginning of cm_rep_cases (to start at 0)
    for i in range(len(cm_rep_cases_dates)-1):
        dates.append(cm_rep_cases_dates[i]-delta)

    cumulative_reported_cases = Bunch(dates=cm_rep_cases_dates, values=cm_rep_cases)
    SEIIIR = Bunch(dates=dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=T)
    return Bunch(cumulative_reported_cases=cumulative_reported_cases, SEIIIR=SEIIIR)

    
