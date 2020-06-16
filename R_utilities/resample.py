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
          p2 <- mean(window_data)
          p3 <- sd(window_data)
        } else {
          p2 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[2])
          p3 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[1])
        }
        c(p2, p3)
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
        samples_negbin = pd.DataFrame(np.zeros((len(idx_range), n_samples)))
    else:
        initial = dat.Confirmed[0]
        daily_increases = np.array(dat.Confirmed[1:dat.shape[0]] - dat.Confirmed[0:(dat.shape[0] - 1)].values)
        daily = np.concatenate(([initial], daily_increases))
        samples_negbin = pd.DataFrame(columns=['s' + str(i) for i in range(1, 101)])
        r = 0
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
            params = r_fit_negbin(window_data, low)
        samples_negbin.loc[r] = MASS.rnegbin(n_samples, params[0], params[1])
        r += 1
    return samples_negbin


# for (j in 1:n_samples){
#     new_df < - dat[range,]
# new_df$Confirmed < - cumsum(samples_negbin[, j])
# # Insert one day of zeroes before so that Bill's code works for all counties
# new_df < - rbind(dat[window,], new_df)
# new_df$Confirmed[1] < - 0
# new_df$FIPS < - as.character(new_df$FIPS)
# new_df$FIPS < - str_pad(new_df$FIPS, 5, pad = '0')
#
# # Where to write these files
# write_path < - paste0('./formatted_data/county_data_resample/', folder_name, '/sample', str_pad(as.character(
#     j), 3, pad = '0'))
# if (!dir.exists(file.path(write_path))){
# dir.create(file.path(write_path))
# }
# write.csv(new_df, file=paste0(write_path, '/', countyfile), row.names = F)


# test = sample_county_negbin(county_files[1051])
# print(test.tail())

for i in county_files:
    resample = sample_county_negbin(i)
    # This returns a dataframe - what do we want to do with it?