import pandas as pd
import numpy as np
from datetime import datetime, timedelta
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
#import math
from pyutilib.misc import timing
from .stochastic import populate_compartments

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri


#main_dir = '../../covid-data/formatted_data/county_data/'
#folder_name = sorted(os.listdir(main_dir))[-1]
#path = main_dir + '/' + folder_name
#
#write_path = '../../covid-data/formatted_data/county_data_resample/' + folder_name + '/'
#if not os.path.exists(write_path):
#    os.makedirs(write_path)
#county_files = sorted(os.listdir(path))

rpy2.robjects.numpy2ri.activate()
rstring="""
    function(window_data, low){
        library(MASS)
        fit <- TRUE
        fit <- tryCatch(fitdistr(window_data, 'Negative Binomial', lower = low), 
                        error = function(cond){
                          return(fitted = FALSE)
                        })
        if (length(fit) == 1){
          p1 <- mean(window_data)
          p2 <- sd(window_data)
        } else {
          p1 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[2])
          p2 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[1])
        }
        c(p1, p2)
    }
"""
r_fit_negbin = robjects.r(rstring)

# Symmetric window, so this is the number of days on either side of the day being calculated,
# leaving a total window size of 2xwindow + 1 days
def sample_county_negbin(*, dat, county, parameters, window=3, n_samples=1):
    # Q: just do this once?
    base = importr('base')
    utils = importr('utils')
    try:
        MASS = importr('MASS')
    except rpy2.robjects.packages.PackageNotInstalledError:
        utils.install_packages('MASS', type='source')
        MASS = importr('MASS')
        

    idx_range = list(range(window, dat.shape[0]))  # Instead using a symmetric window until the end, then using past data

    # If the county has no cases, keep them all at zero
    if dat[county].iloc[-1] == 0:
        return [0] * len(idx_range)
        #samples_negbin = pd.DataFrame(np.zeros((len(idx_range), n_samples)))
    else:
        initial = dat[county][0]
        daily_increases = np.array(dat[county][1:dat.shape[0]] - dat[county][0:(dat.shape[0] - 1)].values)
        daily = np.concatenate(([initial], daily_increases))
        samples_negbin = pd.DataFrame(columns=['s' + str(i) for i in range(1, n_samples+1)])
        #
        # TODO - review this hack.  The reconstruction logic requires the first time step to have
        #           no cases.  Is this reasonable, even if we are resampling?
        #
        samples_negbin.loc[0] = 0
        r = 1
        for i in range((window + 1), dat.shape[0]):
            if i > dat.shape[0] - window:
                window_data = daily[(len(daily) - (2 * window + 1)): len(daily)]
            else:
                # Using a symmetric window (window size is number of days on either side of date of interest)
                window_data = daily[(i - window):(i + window)]
            if (all(window_data == 0)):
                # Need to force the negative binomial parameters to get a fit in some cases
                params = [0,1]
            else:
                if min(window_data) < 2:
                    low = 0.1
                else:
                    low = 1
                if county not in parameters:
                    parameters[county] = {}
                if r in parameters[county]:
                    params = parameters[county][r]
                else:
                    params = parameters[county][r] = r_fit_negbin(window_data, low)
            samples_negbin.loc[r] = MASS.rnegbin(n_samples, params[0], params[1])
            r += 1

    # Reformat: add date and FIPS column, and make a cumulative sum instead of daily counts
    samples_cases = samples_negbin.cumsum()
    return samples_cases['s1'].to_list()

    #print("Z",len(list(dat.index)))
    #Date = pd.DataFrame(dat.index[(window + 1):dat.shape[0]]).reset_index()
    #FIPS = pd.DataFrame(dat.FIPS[(window + 1):dat.shape[0]]).reset_index()
    #samples_cases['Date'] = Date.Date
    #samples_cases['FIPS'] = FIPS.FIPS
    #return samples_cases


#for i in county_files:
#    resample = sample_county_negbin(i)
#    # This returns a dataframe - what do we want to do with it?

def resampled_reconstruction(*, dates, reported_cases_per_day, population, n_steps_per_day,
                              reporting_delay_mean=8, reporting_delay_dev=1.35,
                              reporting_multiplier=10,
                              fixed_incubation=5.2,
                              infectious_lower=2.6, infectious_upper=6, seed=0):
    """
    This function computes a reconstruction of the SEIIIR compartments based on the reported cases.
    1. we resampled the reported cases using an estimated inverse binomial and the reporting_multiplier

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

    return populate_compartments(tcases_timestep=tcases_timestep, 
                                output_reported_cases=output_reported_cases,
                                population=population,
                                fixed_incubation=fixed_incubation,
                                infectious_lower=infectious_lower,
                                infectious_upper=infectious_upper,
                                n_steps_per_day=n_steps_per_day,
                                t_daily_dates=t_daily_dates
                                )

