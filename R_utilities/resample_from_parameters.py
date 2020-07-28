import rpy2.robjects as robjects
import os
from rpy2.robjects.packages import importr
import pandas as pd
import numpy as np
# import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate()

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


county_files = sorted(os.listdir(path))
param_path = '../../covid-data/formatted_data/resample_parameters/'
param_file_name = sorted(os.listdir(param_path))[-1]
paramfile = pd.read_csv(param_path + param_file_name)

n_samples = 1

def draw_county_sample(countyfile, paramfile):
    dat = pd.read_csv(path + '/' + countyfile)
    county_params = paramfile[paramfile['FIPS'] == dat.FIPS[0]]
    # x = 3
    for date in county_params.Date:
        param1 = county_params[county_params.Date == date].param1.values[0]
        param2 = county_params[county_params.Date == date].param2.values[0]
        sample_cases = MASS.rnegbin(n_samples, param1, param2)
        # sample_cases = x
        dat.loc[dat.Date == date, 'Confirmed'] = sample_cases
        # x += 1
    dat = dat[dat.Date <= date]
    dat.Confirmed = dat.Confirmed.cumsum()
    return dat

# Quick test
countyfile = county_files[2000]
newfile = draw_county_sample(countyfile, paramfile)
