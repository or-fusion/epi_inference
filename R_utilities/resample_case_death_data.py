import rpy2.robjects as robjects
import os
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

base = importr('base')
utils = importr('utils')
# utils.install_packages('MASS', type='source')
MASS = importr('MASS')

main_dir = '../../covid-data/formatted_data/county_data/'
folder_name = sorted(os.listdir(main_dir))[-1]
path = main_dir + '/' + folder_name

write_path = '../../covid-data/formatted_data/county_data_resample/' + folder_name + '/'
if not os.path.exists(write_path):
    os.makedirs(write_path)

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

county_files = sorted(os.listdir(path))
# Symmetric window, so this is the number of days on either side of the day being calculated,
# leaving a total window size of 2xwindow + 1 days
window = 3
n_samples = 100

def sample_county_negbin(countyfile):
    dat = pd.read_csv(path + '/' + countyfile)
    idx_range = list(range((window + 1), dat.shape[0]))  # Instead using a symmetric window until the end, then using past data

    # If the county has no cases, keep them all at zero
    if dat.Confirmed.iloc[-1] == 0:
        samples_cases = pd.DataFrame(np.zeros((len(idx_range), n_samples)))
    else:
        # Set up data for sampling case counts
        initial = dat.Confirmed[0]
        daily_increases = np.array(dat.Confirmed[1:dat.shape[0]] - dat.Confirmed[0:(dat.shape[0] - 1)].values)
        daily = np.concatenate(([initial], daily_increases))
        samples_cases = pd.DataFrame(columns=['s' + str(i) for i in range(1, 101)])
        # Set up data for sampling death counts
        initial = dat.Deaths[0]
        daily_increases_death = np.array(dat.Deaths[1:dat.shape[0]] - dat.Deaths[0:(dat.shape[0] - 1)].values)
        daily_deaths = np.concatenate(([initial], daily_increases_death))
        samples_deaths = pd.DataFrame(columns=['s' + str(i) for i in range(1, 101)])
        r = 0
    for i in range((window + 1), dat.shape[0]):
        if i > dat.shape[0] - window:
            window_data = daily[(len(daily) - (2 * window + 1)): len(daily)]
            window_data_deaths = daily_deaths[(len(daily_deaths) - (2 * window + 1)): len(daily_deaths)]
        else:
            # Using a symmetric window (window size is number of days on either side of date of interest)
            window_data = daily[(i - window):(i + window)]
            window_data_deaths = daily_deaths[(i - window):(i + window)]
        if all(window_data == 0):
            # Need to force the negative binomial parameters to get a fit in some cases
            params = [0,1]
        else:
            if min(window_data) < 2:
                low = 0.1
            else:
                low = 1
            params = r_fit_negbin(window_data, low)
        samples_cases.loc[r] = MASS.rnegbin(n_samples, params[0], params[1])

        # Now sample deaths as well:
        if dat.Deaths.iloc[-1] == 0:
            samples_deaths.loc[r] = np.zeros(n_samples)
        else:
            if all(window_data_deaths == 0):
                params = [0,1]
            else:
                if min(window_data_deaths) < 2:
                    low = 0.1
                else:
                    low = 1
                params = r_fit_negbin(window_data_deaths, low)
            samples_deaths.loc[r] = MASS.rnegbin(n_samples, params[0], params[1])
        r += 1

    # Reformat: add date and FIPS column, and make a cumulative sum instead of daily counts
    samples_cases = samples_cases.cumsum()
    samples_deaths = samples_deaths.cumsum()
    Date = pd.DataFrame(dat.Date[(window + 1):dat.shape[0]]).reset_index()
    FIPS = pd.DataFrame(dat.FIPS[(window + 1):dat.shape[0]]).reset_index()
    samples_cases['Date'] = Date.Date
    samples_cases['FIPS'] = FIPS.FIPS
    samples_deaths['Date'] = Date.Date
    samples_deaths['FIPS'] = FIPS.FIPS
    return samples_cases, samples_deaths


cases, deaths = sample_county_negbin(county_files[2040])
print(cases.shape)
print(cases.tail())
print('-------')
print(deaths.shape)
print(deaths.tail())

# for i in county_files:
#     sampled_cases, sampled_deaths = sample_county_negbin(i)
    # This returns two dataframes - what do we want to do with them?