from datetime import datetime, timedelta
import epi_inference.reconstruction.common as common
import numpy as np
from pyutilib.misc.misc import Bunch

"""
This module provides methods for simple reconstruction of the populations
in the compartments from cumulative reported cases. Please see the 
documentation for the individual methods below.
"""

def _transmissions_from_reported_cases(*, dates, reported_cases_per_day, reporting_factor, report_delay):
    """
    This function computes the list of transmissions (new cases) from the reported
    cases, reporting factor, and the reporting delay. The delay doesn't impact the
    number of transmissions, but does impact the estimated dates those transmissions
    occur.

    Parameters
    ----------
    dates : list of python datetime objects
       The days corresponding to the reported cases as datetime objects
    reported_cases_per_day : list of numbers
       The number of reported cases within each day
    reporting_factor : number
       This is the reporting factor. E.g., if reporting_factor=8, then 1 in 8 
       cases are reported
    report_delay : int
       This is the number of days between infection and reporting
    """
    assert len(dates) == len(reported_cases_per_day)
    transmissions = [reporting_factor * r for r in reported_cases_per_day]
    report_delta = timedelta(days=report_delay)
    tdates = list()
    for dt in dates:
        tdates.append(dt - report_delta)
    return Bunch(dates=tdates, values=transmissions)

def _np_transmissions_from_reported_cases(*, dates, counties, reported_cases_per_day, reporting_factor, report_delay):
    """
    This function computes the list of transmissions (new cases) from the reported
    cases, reporting factor, and the reporting delay. The delay doesn't impact the
    number of transmissions, but does impact the estimated dates those transmissions
    occur.

    Parameters
    ----------
    dates : numpy array of datetime objects
        array of datetime objects corresponding to the rows of reported_cases_per_day
    counties : numpy chararray
       the names of the counties (or nodes) corresponding to the coumns in the
       reported_cases_per_day
    reported_cases_per_day : Numpy two-dimensional array
       This is a numpy array that contains the reported cases per day. Each row
       corresponds to a different day, and each column is a different county (node).
    reporting_factor : number
       This is the reporting factor. E.g., if reporting_factor=8, then 1 in 8 
       cases are reported
    report_delay : int
       This is the number of days between infection and reporting
    """
    # check the types
    assert isinstance(dates, np.ndarray) and dates.dtype == np.object
    assert isinstance(counties, np.ndarray) and counties.dtype == np.object
    
    assert len(dates) == reported_cases_per_day.shape[0]
    assert len(counties) == reported_cases_per_day.shape[1]
    transmissions = reporting_factor * reported_cases_per_day
    report_delta = timedelta(days=report_delay)
    transmission_dates = np.asarray([dates[i] - report_delta for i in range(len(dates))])
    return Bunch(dates=transmission_dates, values=transmissions)


def reconstruct_states_deterministic_decay(*, dates, reported_cases_per_day, population, sigma, gamma, reporting_factor, report_delay, county=None, warnings=None):
    """
    This function reconstructs the state of the system using the
    cumulative reported cases and assuming a constant reporting
    delay. Using this and the reporting factor, we can compute the
    approximate transmissions each day.  These are then used to
    reconstruct the other states, assuming a discrete time model with
    a simple decay term for leaving the E, I1, I2, and I3
    compartments. The model is discretized by days.
    
    Parameters
    ----------
    dates : list
        list of datetime objects corresponding to the dates of the reported_cases_per_day

    reported_cases_per_day : list
        list of the reported cases per day

    population : number
       This is the overall population

    sigma : float
       This is the rate of transfer out of the exposed compartment (e.g., 1/(incubation period))

    gamma : float
       This model includes 3 infectious compartments. Therefore this is an approximation of the 
       overall rate of transfer out of the I compartment. Here, we create three compartments,
       each with a rate of 3*gamma (where gamma is 1/(infectious period))

    reporting_factor : float
       This factor accounts for under-reporting. If reporting_factor=8, then 1 in 8 cases are reported

    report_delay : int
       This is the number of days between infection and reporting

    county : string
       The county name, used for debugging output

    warnings : list
       A list that can be used to store warnings for the user

    Returns
    -------
       tuple : (dates, T, S, E, I1, I2, I3, R)
          dates: dates corresponding to the states in the model
          T: list of new transmissions
          S: list of susceptible population values
          E: list of exposed population values
          I1, I2, I3: lists of infective counts in I1, I2, and I3 compartments respectively
          R: list of recovered population values
    """
    assert population > 1
    assert sigma >=0
    assert gamma >= 0
    assert reporting_factor >= 1
    assert report_delay >= 1 and type(report_delay) is int

    transmissions = _transmissions_from_reported_cases(dates=dates,
                                            reported_cases_per_day=reported_cases_per_day,
                                            reporting_factor=reporting_factor,
                                            report_delay=report_delay)
    transmission_dates = transmissions.dates
    transmissions = transmissions.values

    # create lists to store compartment numbers
    S = [None]*len(transmissions)
    E = [None]*len(transmissions)
    I1 = [None]*len(transmissions)
    I2 = [None]*len(transmissions)
    I3 = [None]*len(transmissions)
    R = [None]*len(transmissions)

    # assume fully susceptible population to start
    S[0] = population
    E[0] = 0
    I1[0] = 0
    I2[0] = 0
    I3[0] = 0
    R[0] = 0

    # main simulation loop
    for t in range(len(transmissions)-1):
        delta0 = transmissions[t]
        if delta0 > S[t] and not warnings is None:
            warnings.append("WARNING: Cases in county %s exceeded population size (%f > %f) at time step %d. Ignoring cases that exceed population size." % (str(county), delta0, S[t], t))
            delta0 = S[t]
        S[t+1] = S[t] - delta0
        #if S[t+1] <= 0:
        #    # print(sum(reported_cases_per_day))
        #    # print(cumulative_reported_cases[-1])
        #    # print(population)
        #    # for idx,d in enumerate(tdates):
        #    #      print(d,':',S[idx], '->', transmissions[idx])
        #    raise ValueError("reconstruct_states_deterministic_decay computed a negative"
        #                     " value for the susceptible population. This likely means that"
        #                     " the reported cases or the reporting_factor are too large"
        #                     " for the population value specified. This happened at "
        #                     "timestep: {} with date: {}, transmissions: {}, and S[t+1]: {}".format(
        #                         t, rdates[t], transmissions[t], S[t+1])
        #                     )
        
        delta1 = min(sigma*E[t], E[t] + delta0)
        E[t+1] = E[t] + delta0 - delta1
        #E[t+1] = E[t] + transmissions[t] - sigma*E[t]
        #E[t+1] = max(0, E[t+1])
        assert E[t+1] >= 0

        delta2 = min(gamma*3*I1[t], I1[t] + delta1)
        I1[t+1] = I1[t] + delta1 - delta2
        #I1[t+1] = I1[t] + sigma*E[t] - gamma*3*I1[t]
        #I1[t+1] = max(0, I1[t+1])
        assert I1[t+1] >= 0
        
        delta3 = min(gamma*3*I2[t], I2[t] + delta2)
        I2[t+1] = I2[t] + delta2 - delta3
        #I2[t+1] = I2[t] + gamma*3*I1[t] - gamma*3*I2[t]
        #I2[t+1] = max(0, I2[t+1])
        assert I2[t+1] >= 0
        
        delta4 = min(gamma*3*I3[t], I3[t] + delta3)
        I3[t+1] = I3[t] + delta3 - delta4
        #I3[t+1] = I3[t] + gamma*3*I2[t] - gamma*3*I3[t]
        #I3[t+1] = max(0, I3[t+1])
        assert I3[t+1] >= 0
        
        R[t+1] = R[t] + delta4
        #R[t+1] = R[t] + gamma*3*I3[t]
        assert R[t+1] >= 0

    return Bunch(dates=transmission_dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=transmissions)

def np_reconstruct_states_deterministic_decay(*, dates, counties, reported_cases_per_day, populations, sigma, gamma, reporting_factor, report_delay):
    """
    This function reconstructs the state of the system using the
    cumulative reported cases and assuming a constant reporting
    delay. Using this and the reporting factor, we can compute the
    approximate transmissions each day.  These are then used to
    reconstruct the other states, assuming a discrete time model with
    a simple decay term for leaving the E, I1, I2, and I3
    compartments. The model is discretized by days.
    
    Parameters
    ----------
    dates : numpy array of datetime objects
        array of datetime objects corresponding to the rows of reported_cases_per_day
    counties : numpy array of objects (strings)
       the names of the counties (or nodes) corresponding to the coumns in the
       reported_cases_per_day
    reported_cases_per_day : Numpy two-dimensional array
       This is a numpy array that contains the reported cases per day. Each row
       corresponds to a different day, and each column is a different county (node).
    populations : numpy array of populations
       This is an array of populations. Each entry corresponds to one of the counties
    sigma : float
       This is the rate of transfer out of the exposed compartment (e.g., 1/(incubation period))
    gamma : float
       This model includes 3 infectious compartments. Therefore this is an approximation of the 
       overall rate of transfer out of the I compartment. Here, we create three compartments,
       each with a rate of 3*gamma (where gamma is 1/(infectious period))
    reporting_factor : float
       This factor accounts for under-reporting. If reporting_factor=8, then 1 in 8 cases are reported

    report_delay : int
       This is the number of days between infection and reporting

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
    
    ndates = len(dates)
    ncounties = len(counties)
    assert reported_cases_per_day.shape[0] == ndates
    assert reported_cases_per_day.shape[1] == ncounties
    
    assert np.all(populations > 0)
    assert sigma >=0
    assert gamma >= 0
    assert reporting_factor >= 1
    assert report_delay >= 1 and type(report_delay) is int

    transmissions = \
        _np_transmissions_from_reported_cases(dates=dates,
                                              counties=counties,
                                              reported_cases_per_day=reported_cases_per_day,
                                              reporting_factor=reporting_factor,
                                              report_delay=report_delay)
    transmission_dates = transmissions.dates
    transmissions = transmissions.values

    # create arrays to store compartment numbers
    S  = np.NaN*np.zeros((ndates,len(counties)))
    E  = np.NaN*np.zeros((ndates,len(counties)))
    I1 = np.NaN*np.zeros((ndates,len(counties)))
    I2 = np.NaN*np.zeros((ndates,len(counties)))
    I3 = np.NaN*np.zeros((ndates,len(counties)))
    R  = np.NaN*np.zeros((ndates,len(counties)))

    # assume fully susceptible population to start
    S[0,:] = populations
    E[0,:] = 0
    I1[0,:] = 0
    I2[0,:] = 0
    I3[0,:] = 0
    R[0,:] = 0
    
    # main simulation loop
    for t in range(ndates-1):
        delta0 = transmissions[t,:]
        if np.any(delta0 > S[t,:]):
            print("WARNING: cases exceeded population size.  Ignoring cases that exceed population size.")
            delta0 = np.minimum(delta0, S[t,:])
        S[t+1,:] = S[t,:] - delta0
        
        delta1 = np.minimum(sigma*E[t,:], E[t,:] + delta0)
        E[t+1,:] = E[t,:] + delta0 - delta1

        delta2 = np.minimum(gamma*3*I1[t,:], I1[t,:] + delta1)
        I1[t+1,:] = I1[t,:] + delta1 - delta2
        
        delta3 = np.minimum(gamma*3*I2[t,:], I2[t,:] + delta2)
        I2[t+1,:] = I2[t,:] + delta2 - delta3
        
        delta4 = np.minimum(gamma*3*I3[t,:], I3[t,:] + delta3)
        I3[t+1,:] = I3[t,:] + delta3 - delta4
        
        R[t+1,:] = R[t,:] + delta4

    assert np.any(np.isfinite(S))
    assert np.any(np.isfinite(E))
    assert np.any(np.isfinite(I1))
    assert np.any(np.isfinite(I2))
    assert np.any(np.isfinite(I3))
    assert np.any(np.isfinite(R))
    assert np.any(np.isfinite(transmissions))
    
    return Bunch(dates=transmission_dates, S=S, E=E, I1=I1, I2=I2, I3=I3, R=R, transmissions=transmissions)

