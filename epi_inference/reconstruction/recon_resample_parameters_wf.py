import sys
import os
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

from ..engine.task import Task
from ..engine.task_registry import register_task


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

def initialize_R():
    rpy2.robjects.numpy2ri.activate()
    base = importr('base')
    utils = importr('utils')
    try:
        MASS = importr('MASS')
    except:
        utils.install_packages('MASS', type='source')
        MASS = importr('MASS')
    r_fit_negbin = robjects.r(rstring)


# Symmetric window, so this is the number of days on either side of the day being calculated,
# leaving a total window size of 2xwindow + 1 days
window = 3
# n_samples = 100

def save_county_negbin_parameters_cases(data_directory):
    county_files = sorted(os.listdir(data_directory))
    param_cases = pd.DataFrame()
    # Set up files to store the data
    for countyfile in county_files:
        county_case_params = pd.DataFrame()
        dat = pd.read_csv(data_directory + '/' + countyfile, encoding='latin-1')
        county_case_params['FIPS'] = pd.Series(dat.FIPS[(window + 1):dat.shape[0]].values)
        county_case_params['Date'] = pd.Series(dat.Date[(window + 1):dat.shape[0]].values)
        county_case_params['param1'] = pd.Series([0] * county_case_params.shape[0])
        county_case_params['param2'] = pd.Series([1] * county_case_params.shape[0])
        idx_range = list(range((window + 1), dat.shape[0]))  # Instead using a symmetric window until the end, then using past data
        # If the county has no cases, keep them all at zero
        if dat.Confirmed.iloc[-1] == 0:
            # We're done - add county_case_params to param_cases dataframe and move on.
            param_cases = param_cases.append(county_case_params)
        else:
            # Set up data for parameter fitting
            initial = dat.Confirmed[0]
            daily_increases = np.array(dat.Confirmed[1:dat.shape[0]] - dat.Confirmed[0:(dat.shape[0] - 1)].values)
            daily = np.concatenate(([initial], daily_increases))
            r = 0
            for i in range((window + 1), dat.shape[0]):
                if i > dat.shape[0] - window:
                    window_data = daily[(len(daily) - (2 * window + 1)): len(daily)]
                else:
                    # Using a symmetric window (window size is number of days on either side of date of interest)
                    window_data = daily[(i - window):(i + window)]
                if all(window_data == 0):
                    # Need to force the negative binomial parameters to get a fit in some cases
                    county_case_params.loc[r, 'param1'] = 0
                    county_case_params.loc[r, 'param2'] = 1
                else:
                    if min(window_data) < 2:
                        low = 0.1
                    else:
                        low = 1
                    params = r_fit_negbin(window_data, low)
                    county_case_params.loc[r, 'param1'] = params[0]
                    county_case_params.loc[r, 'param2'] = params[1]
                r += 1
            param_cases = param_cases.append(county_case_params)
    return param_cases


class ReconstructionEstimateParameters(Task):

    def __init__(self):
        Task.__init__(self, "reconstruction_estimate_parameters",
            "Estimate parameters used for resampling.")

    def validate(self, CONFIG):
        valid_options = set(['data_directory', 'output_csv', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        initilize_R()
        self._warnings = []
        self.validate(CONFIG)
        run(CONFIG, self._warnings)
        save_params = save_county_negbin_parameters_cases(CONFIG['data_directory'])
        save_params.to_csv(CONFIG['output_csv'], index=False)


register_task(ReconstructionEstimateParameters())


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('You must provide the path to the data directory')
        quit()

    fname = sys.argv[1]
    print(fname)
    save_params = save_county_negbin_parameters_cases(fname)
    save_params.to_csv('../../covid-data/formatted_data/resample_parameters/negbin_params_' + fname[-11:-1] + '.csv', index=False)
