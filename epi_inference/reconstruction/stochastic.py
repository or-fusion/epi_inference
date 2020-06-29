import pandas as pd
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
from pyutilib.misc.misc import Bunch

def stochastic_reconstruction(*, dates, reported_cases_per_day, population, n_steps_per_day,
                              reporting_delay_mean=8, reporting_delay_dev=1.35,
                              reporting_multiplier=10,
                              fixed_incubation=5.2,
                              infectious_lower=2.6, infectious_upper=6, seed=0):
    """
    This function computes a reconstruction of the SEIIIR compartments based on the reported cases.
    1. for each day, the actual number of cases is sampled from an inverse binomial using the reported
        cases and the reporting_multiplier
    2. for each case, we sample the reporting_delay using a log-normal with the reporting_delay_mean and
       reporting_delay_dev - this produces the vector of transmissions in time
    3. The transmissions in time cause movement from S->E, and then a discrete-time SEIIIR model is used
       to populate the compartments.

    Parameters
    ----------
    dates : list
        list of datetime objects corresponding to the dates of the reported_cases_per_day
    reported_cases_per_day : list
        list of the reported cases per day
    population : float
        population of this node (county)
    n_steps_per_day : int
        the number of timesteps to take in one day (e.g., value of 4 means the delta_t = 0.25 days)
    reporting_delay_mean : float
        reporting delay is drawn from a log normal
           d=np.random.lognormal(np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev))
        Note: If reporting_delay_dev is set to None, then the reporting delay is fixed to the value of
        reporting_delay_mean and not sampled.
    reporting_delay_dev : float
        reporting delay is drawn from a log normal
           d=np.random.lognormal(np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev))
        Note: If this is set to None, then the reporting delay is fixed to the value of reporting_delay_mean
        and not sampled.
    reporting_multiplier : int
        the number of actual cases per reported case (i.e., reporting fraction is 1/reporting_multiplier)
    fixed_incubation : float
        the incubation period in days
    infectious_lower : float
        the infectious period is drawn from a uniform distribution between infectious_lower and infectious_upper
    infectious_upper : float
        the infectious period is drawn from a uniform distribution between infectious_lower and infectious_upper
    seed: int
        the seed for the random number generator

    Returns
    -------
        dict : a dict of the dates and the compartments
    """
    if seed:
        np.random.seed(seed)

    n_r_days = len(dates)

    # create a vector of the infections in time (when the infection occurred - S->E)
    # this is called T (transmissions) in the code -> see tdates, tcases
    # this needs to be padded at the front with days prior to the reporting dates
    # and converted to timesteps that includes "steps per day"
    assert type(n_steps_per_day) is int
    assert n_steps_per_day >= 1
    padding_days = 50 # days - maybe this should be passed in?
    padding_timesteps = padding_days*n_steps_per_day
    tcases_timestep = [0]*(n_r_days+padding_days)*n_steps_per_day

    # probability of confirmation of case - reporting fraction
    # this should be drawn from a distribution as well if good data exists
    p = 1 / reporting_multiplier

    # loop through each day and draw the total number of cases that day
    # from the reported cases and the reporting fraction
    # These are the total number of cases that could be reportable
    for r_day in range(len(dates)):
        r_timestep = r_day*n_steps_per_day
        if reported_cases_per_day[r_day] > 0:
            # draw the total number of reportable cases
            reportable_cases_day = int(reported_cases_per_day[r_day] + np.random.negative_binomial(reported_cases_per_day[r_day],p))
#        else:
#            # draw the total number of reportable cases
#            reportable_cases_day = max(0,np.random.negative_binomial(1,p)-1)
#            print('Zero reported cases on date {} ({}). Reportable cases = {}'.format(dates[r_day], r_day, reportable_cases_day))
#
#        if True:
            if reportable_cases_day > 0:
                # now draw the delays from infection to confirmation (log normal)
                # one delay is drawn for each reportable case and is in units of days
                if reporting_delay_dev is None:
                    # used for testing against simulated data
                    delays_days = [reporting_delay_mean]*reportable_cases_day
                    delays_days = [d*n_steps_per_day for d in delays_days]
                else:
                    delays_days = np.random.lognormal(mean=np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev), size=reportable_cases_day)
                    # convert this to timesteps
                    assert type(delays_days) is np.ndarray
                    delays_days = delays_days*n_steps_per_day
                delays_timesteps = np.round(delays_days).astype(int)
                # add one transmission (infection) based on each reportable case and its delay
                for delay_timestep in delays_timesteps:
                    if r_timestep - delay_timestep + padding_timesteps >= 0:
                        tcases_timestep[r_timestep - delay_timestep + padding_timesteps] += 1

    # truncate 2*reporting_delay_mean days off of the transmissions
    # since we don't yet have the reported cases to estimate this appropriately
    # 2*reporting_delay_mean might not be enough
    int_delay = int(np.round(2*reporting_delay_mean))
    t_daily_dates = [dates[0] + timedelta(days=i) for i in range(-padding_days, n_r_days - int_delay)]
    tcases_timestep = tcases_timestep[: -int_delay * n_steps_per_day]

    # capture the reported cases per day - add the padding and truncate the end
    output_reported_cases = [0]*padding_days
    output_reported_cases.extend(reported_cases_per_day[: -int_delay])

    # use these transmissions to generate a single stochastic reconstruction
    # of each of the compartments (SEIIIR)
    S = np.zeros(len(tcases_timestep))
    S[0] = population
    E = np.zeros(len(tcases_timestep))
    I1 = np.zeros(len(tcases_timestep))
    I2 = np.zeros(len(tcases_timestep))
    I3 = np.zeros(len(tcases_timestep))
    R = np.zeros(len(tcases_timestep))

    # Approach 1:
    # - sigma is fixed at 1/5.2 days
    # - serial interval is drawn from uniform 6.5-8.2, and used
    #   to compute gamma
    # It might be better to just draw sigma and gamma from separate
    # distributions here.
    # Here, the serial interval is drawn once and that value is used
    # for the entire simulation, however, it could be drawn for each day
    # as well.
    # Note: These are drawn for every day, but they are the same for each day
    # - this is done so the code is ready for different values on each day if
    #   we decide to do that.
    sigma = 1.0/fixed_incubation*np.ones(len(tcases_timestep)) # days
    infectious_period = np.random.uniform(infectious_lower, infectious_upper)*np.ones(len(tcases_timestep)) # CDL , size=len(tcases_timestep))
    gamma = np.reciprocal(infectious_period)
    #incubation = np.random.exponential(scale=mean_incubation, size=len(tcases))
    # CDL serial_interval = np.random.uniform(6.35, 8.05)*np.ones(len(tcases_timestep))
    #serial_interval = np.random.uniform(6.5, 8.2)*np.ones(len(tcases_timestep))
    # compute gamma = 1/ ( 2*(SI-1/sigma) )
    #temp = np.copy(sigma)
    #np.reciprocal(temp, out=temp)
    #temp = 2*(serial_interval-temp)
    #gamma = np.reciprocal(temp)
    
    prob_E = 1-np.exp(-1.0/n_steps_per_day*sigma)
    prob_III = 1-np.exp(-1.0/n_steps_per_day*3*gamma)

    # loop through all of the days and compute the compartments
    # with the stochastic simulations
    for t in range(len(tcases_timestep)-1):
        Stout = tcases_timestep[t]
        if Stout >= S[t]:
            # reconstruction indicates the number of infections
            # exceed the population - flag this for error reporting
            # means the reported cases or reporting factor are too large
            # for the population value specified
            # Todo: flag a warning
            Stout = S[t]
            
        S[t+1] = S[t] - Stout
        Etout = np.random.binomial(E[t], prob_E[t])
        E[t+1] = E[t] + Stout - Etout
        I1tout = np.random.binomial(I1[t], prob_III[t])
        I1[t+1] = I1[t] + Etout - I1tout
        I2tout = np.random.binomial(I2[t], prob_III[t])
        I2[t+1] = I2[t] + I1tout - I2tout
        I3tout = np.random.binomial(I3[t], prob_III[t])
        I3[t+1] = I3[t] + I2tout - I3tout
        R[t+1] = R[t] + I3tout

    # now bring these vectors back to days
    S = [S[t] for t in range(len(S)) if t % n_steps_per_day == 0]
    E = [E[t] for t in range(len(E)) if t % n_steps_per_day == 0]
    I1 = [I1[t] for t in range(len(I1)) if t % n_steps_per_day == 0]
    I2 = [I2[t] for t in range(len(I2)) if t % n_steps_per_day == 0]
    I3 = [I3[t] for t in range(len(I3)) if t % n_steps_per_day == 0]
    R = [R[t] for t in range(len(R)) if t % n_steps_per_day == 0]
    T = [None]*len(S)
    for day_idx in range(len(T)):
        T[day_idx] = 0
        for t in range(n_steps_per_day):
            T[day_idx] += tcases_timestep[t+day_idx*n_steps_per_day]
    assert len(t_daily_dates) == len(T)
    assert len(t_daily_dates) == len(S)
    # return the output
    # - rdates, rcases: reported cases and their dates
    #   (new reported cases each day - not cumulative)
    # - tdates: dates for all the other compartments
    # - tcases: transmissions (reportable cases at the time
    #   of the initial transmission
    # - S,E,I1,I2,I3,R: numbers of individuals in each of the
    #   compartments on the dates in tdates
    return Bunch(dates=t_daily_dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=T, orig_rep_cases=output_reported_cases)

def np_stochastic_reconstruction(*, dates,
                                 reported_cases_per_day,
                                 counties,
                                 populations,
                                 n_steps_per_day,
                                 reporting_delay_mean=8, reporting_delay_dev=1.35,
                                 reporting_multiplier=10,
                                 fixed_incubation=5.2,
                                 infectious_lower=2.6, infectious_upper=6):
    """
    This function computes a reconstruction of the SEIIIR compartments based on the reported cases.
    1. for each day, the actual number of cases is sampled from an inverse binomial using the reported
        cases and the reporting_multiplier
    2. for each case, we sample the reporting_delay using a log-normal with the reporting_delay_mean and
       reporting_delay_dev - this produces the vector of transmissions in time
    3. The transmissions in time cause movement from S->E, and then a discrete-time SEIIIR model is used
       to populate the compartments.

    Parameters
    ----------
    dates : numpy array of datetime objects
       The dates corresponding to the rows in the cumulative reported cases
    counties : numpy array of object (strings)
       The names of the counties (or nodes) corresponding to the columns in the
       cumulative reported cases
    reported_cases_per_day : Numpy two-dimensional array
       This is a numpy array that contains the reported cases per day. Each row
       corresponds to a different day, and each column is a different county (node).
    populations : numpy array of populations
       This is an array of populations. Each entry corresponds to one of the counties
    n_steps_per_day : int
        the number of timesteps to take in one day (e.g., value of 4 means the delta_t = 0.25 days)
    reporting_delay_mean : float
        reporting delay is drawn from a log normal
           d=np.random.lognormal(np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev))
        Note: If reporting_delay_dev is set to None, then the reporting delay is fixed to the value of
        reporting_delay_mean and not sampled.
    reporting_delay_dev : float
        reporting delay is drawn from a log normal
           d=np.random.lognormal(np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev))
        Note: If this is set to None, then the reporting delay is fixed to the value of reporting_delay_mean
        and not sampled.
    reporting_multiplier : int
        the number of actual cases per reported case (i.e., reporting fraction is 1/reporting_multiplier)
    fixed_incubation : float
        the incubation period in days
    infectious_lower : float
        the infectious period is drawn from a uniform distribution between infectious_lower and infectious_upper
    infectious_upper : float
        the infectious period is drawn from a uniform distribution between infectious_lower and infectious_upper

    Returns
    -------
       Bunch : (like a dict with keys: dates, S, E, I1, I2, I3, R, transmissions)
          dates: dates corresponding to the states in the model
          S: numpy array of susceptible population values
          E: numpy array of exposed population values
          I1, I2, I3: numpy arrays of infective counts in I1, I2, and I3 compartments respectively
          R: numpy array of recovered population values
          transmissions: numpy array of transmissions (from S->E) 
    """
    # check the types
    assert isinstance(dates, np.ndarray) and dates.dtype == np.object
    assert isinstance(counties, np.ndarray) and counties.dtype == np.object
    assert isinstance(reported_cases_per_day, np.ndarray) and reported_cases_per_day.dtype == np.float
    
    n_r_days = len(dates)
    n_counties = len(counties)
    assert reported_cases_per_day.shape[0] == n_r_days
    assert reported_cases_per_day.shape[1] == n_counties
    assert np.all(populations > 0)

    # create an array of the infections in time (when the infection occurred - S->E)
    # this is called tcases (transmissions) in the code
    # the tcases array needs to be padded at the front with days prior to the reporting dates
    # and converted to timesteps that includes "steps per day"
    assert type(n_steps_per_day) is int
    assert n_steps_per_day >= 1
    padding_days = 30 # days - maybe this should be passed in?
    padding_timesteps = padding_days*n_steps_per_day
    n_r_daily_dates = len(dates) # number of days of reported cases
    n_timesteps = (n_r_days+padding_days)*n_steps_per_day
    tcases_timestep = np.zeros((n_timesteps,n_counties))

    # probability of confirmation of case - reporting fraction
    # this should be drawn from a distribution as well if good data exists
    p = 1 / reporting_multiplier

    # loop through each day and draw the total number of cases that day
    # from the reported cases and the reporting fraction
    # These are the total number of cases that could be reportable
    for r_day in range(len(dates)):
        r_timestep = r_day*n_steps_per_day
        reported_cases = reported_cases_per_day[r_day,:]
        #foo = np.random.negative_binomial([100,10000], p=0.1, size=(1000,2))
        # loop over each of the counties
        for c, cname in enumerate(counties):
            if reported_cases[c] > 0:
                # draw the total number of reportable cases
                reportable_cases_day = int(reported_cases[c] + np.random.negative_binomial(reported_cases[c],p))
                if reportable_cases_day > 0:
                    # now draw the delays from infection to confirmation (log normal)
                    # one delay is drawn for each reportable case and is in units of days
                    if reporting_delay_dev is None:
                        # used for testing against simulated data
                        delays_days = reporting_delay_mean*np.ones(reportable_cases_day)
                        # convert to timesteps instead of days
                        delays_days = n_steps_per_day*delay_days
                    else:
                        delays_days = np.random.lognormal(mean=np.log(reporting_delay_mean), sigma=np.log(reporting_delay_dev), size=reportable_cases_day)
                        # convert this to timesteps instead of days
                        assert type(delays_days) is np.ndarray
                        delays_timesteps = delays_days*n_steps_per_day
                    # round to integer timestep
                    delays_timesteps = np.round(delays_timesteps).astype(int)
                    # add one transmission (infection) based on each reportable case and its delay
                    for delay_timestep in delays_timesteps:
                        tcases_timestep[r_timestep - delay_timestep + padding_timesteps,c] += 1

    # truncate 2*reporting_delay_mean days off of the transmissions
    # since we don't yet have the reported cases to estimate this appropriately
    # 2*reporting_delay_mean might not be enough
    t_daily_dates = np.asarray([dates[0] + timedelta(days=i) for i in range(-padding_days, n_r_days-2*reporting_delay_mean)],dtype=object)
    tcases_timestep = tcases_timestep[:-2*reporting_delay_mean*n_steps_per_day,:]
    n_timesteps = len(tcases_timestep)
    
    # create arrays to store compartment numbers
    S  = np.zeros((n_timesteps,len(counties)), dtype=int)
    E  = np.zeros((n_timesteps,len(counties)), dtype=int)
    I1 = np.zeros((n_timesteps,len(counties)), dtype=int)
    I2 = np.zeros((n_timesteps,len(counties)), dtype=int)
    I3 = np.zeros((n_timesteps,len(counties)), dtype=int)
    R  = np.zeros((n_timesteps,len(counties)), dtype=int)

    # assume fully susceptible population to start
    S[0,:] = populations
    E[0,:] = 0
    I1[0,:] = 0
    I2[0,:] = 0
    I3[0,:] = 0
    R[0,:] = 0

    # Approach 1:
    # - sigma is fixed at 1/5.2 days
    # - serial interval is drawn from uniform 6.5-8.2, and used
    #   to compute gamma
    # It might be better to just draw sigma and gamma from separate
    # distributions here.
    # Here, the serial interval is drawn once and that value is used
    # for the entire simulation, however, it could be drawn for each day
    # as well.
    # Note: These are drawn for every day, but they are the same for each day
    # - this is done so the code is ready for different values on each day if
    #   we decide to do that.
    sigma = 1.0/fixed_incubation # units of days
    infectious_period = np.random.uniform(infectious_lower, infectious_upper)*np.ones(len(tcases_timestep)) # CDL, size=len(tcases_timestep))
    gamma = np.reciprocal(infectious_period)
    #incubation = np.random.exponential(scale=mean_incubation, size=len(tcases))
    # CDL serial_interval = np.random.uniform(6.35, 8.05)*np.ones(len(tcases_timestep))
    #serial_interval = np.random.uniform(6.5, 8.2)*np.ones(len(tcases_timestep))
    # compute gamma = 1/ ( 2*(SI-1/sigma) )
    #temp = np.copy(sigma)
    #np.reciprocal(temp, out=temp)
    #temp = 2*(serial_interval-temp)
    #gamma = np.reciprocal(temp)
    
    prob_E = 1-np.exp(-1.0/n_steps_per_day*sigma)
    prob_III = 1-np.exp(-1.0/n_steps_per_day*3*gamma)

    # loop through all of the days and compute the compartments
    # with the stochastic simulations
    for t in range(len(tcases_timestep)-1):
        Stout = tcases_timestep[t,:]
        if np.any(Stout >= S[t,:]):
            # reconstruction indicates the number of infections
            # exceed the population - flag this for error reporting
            # means the reported cases or reporting factor are too large
            # for the population value specified
            # Todo: flag a warning
            print("WARNING: cases exceeded population size.  Ignoring cases that exceed population size.")
            Stout = np.minimum(Stout, S[t,:])
            
        S[t+1,:] = S[t,:] - Stout
        Etout = np.random.binomial(E[t,:], prob_E, size=n_counties)
        E[t+1,:] = E[t,:] + Stout - Etout
        I1tout = np.random.binomial(I1[t,:], prob_III[t], size=n_counties)
        I1[t+1,:] = I1[t,:] + Etout - I1tout
        I2tout = np.random.binomial(I2[t,:], prob_III[t], size=n_counties)
        I2[t+1,:] = I2[t,:] + I1tout - I2tout
        I3tout = np.random.binomial(I3[t,:], prob_III[t], size=n_counties)
        I3[t+1,:] = I3[t,:] + I2tout - I3tout
        R[t+1,:] = R[t,:] + I3tout

    # now bring these vectors back to days
    timesteps_days = [t for t in range(len(S)) if t % n_steps_per_day == 0] 
    S = S[timesteps_days]
    E = E[timesteps_days]
    I1 = I1[timesteps_days]
    I2 = I2[timesteps_days]
    I3 = I3[timesteps_days]
    R = R[timesteps_days]
    T = np.zeros((len(t_daily_dates)+1, n_counties))
    T[1:] = np.cumsum(tcases_timestep, axis=0)[timesteps_days]
    T = np.diff(T, axis=0)
    assert len(t_daily_dates) == T.shape[0]
    assert n_counties == T.shape[1]
    assert len(t_daily_dates) == S.shape[0]
    assert n_counties == S.shape[1]

    assert np.any(np.isfinite(S))
    assert np.any(np.isfinite(E))
    assert np.any(np.isfinite(I1))
    assert np.any(np.isfinite(I2))
    assert np.any(np.isfinite(I3))
    assert np.any(np.isfinite(R))
    assert np.any(np.isfinite(T))
    
    return Bunch(dates=t_daily_dates, counties=counties, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=T)
