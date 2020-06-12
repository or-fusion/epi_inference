import sys
import pandas as pd
pd.set_option("display.max_rows", None)
import csv
from os import listdir, makedirs
from os.path import isfile, join, splitext, exists
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

for shortlong in ['short', 'long']:
    output_path = './data/florida-{}/'.format(shortlong)
    if not exists(output_path):
        makedirs(output_path)
    for strR0 in ['2.25']: #['1.25', '1.50', '1.75', '2.00', '2.25', '2.50', '2.75', '3.00']:
        path = '../pipeline_simulation/florida/model_output/florida-R0-{}-{}_None'.format(strR0, shortlong)

        # get the list of all the parquet files (one for each realization)
        files = list()
        for f in listdir(path):
            if isfile(join(path, f)):
                filename, extension = splitext(f)
                if extension == '.parquet':
                    files.append(filename)

        # read the parquet files and output the data in csv
        d = dict()
        counties = None
        dates = None
        # loop over all realizations
        for i,f in enumerate(files):
            print('processing realization ... ', i)
            parquet_file = join(path, '{}.parquet'.format(f))
            df = pd.read_parquet(parquet_file, engine='pyarrow')

            # write out the original file
            fname = join(output_path, 'FL_{}_realization_{}_original_parquet.csv'.format(strR0, i))
            df.to_csv(fname, index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")

            # adjust the names of the columns to align with ours
            # change time to  Date
            # change the FIPS codes from float to string
            rename_dict = {'time': 'Date'}
            for c in df.columns:
                if c.endswith('.0'):
                    rename_dict[c] = c[:-2]
            df = df.rename(columns=rename_dict)
            # change the order so Date is first
            reorder = ['Date']
            for c in df.columns:
                if c != 'Date':
                    reorder.append(c)
            df = df.reindex(columns=reorder)

            #  write the  SEIIIR output
            fname = join(output_path, 'FL_SEIIIR_R0_{}_{}_realization_{}.csv'.format(strR0, shortlong, i))
            df.to_csv(fname, index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")

            # build a dataframe of their cumulative reported cases
            dfcmI = df[df['comp'] == 'cumI'].copy()
            del dfcmI['comp']

            # shift all the cases by 3 days (they count cumI at ~5 days and we want a total of 8 days)
            mask = ~(dfcmI.columns.isin(['Date']))
            cols_to_shift = dfcmI.columns[mask]
            dfcmI[cols_to_shift] = dfcmI[cols_to_shift].shift(+3, fill_value=0.0).copy()
            # let's apply a reporting factor of 1/10 and then round
            dfcmI = dfcmI.set_index('Date')
            dfcmI = dfcmI.multiply(0.1)
            dfcmI = dfcmI.round()
            dfcmI = dfcmI.reset_index()
            fname = join(output_path, 'FL_cumulative_I_shifted_R0_{}_{}_realization_{}.csv'.format(strR0, shortlong, i))
            dfcmI.to_csv(fname, index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")

            ###
            # build a dataframe of cumulative reported cases using stochastic
            # process parameters from Derek - get the transmissions from the
            # simulated values in the S compartment
            dftx = df[df['comp'] == 'S']
            del dftx['comp']
            dftx = dftx.set_index('Date')
            dftx = dftx.diff()
            dftx.iloc[0] = 0
            dftx = dftx.abs().astype(int)
            # dftx is now the transmissions with the date that they occured
            # for each of these transmissions, let's find out when they are reported
            dfrepcases = pd.DataFrame({'Date': dftx.index})
            for c in dftx.columns:
                dftxc = dftx[c]
                dftxcv = dftxc.values
                repcases = np.zeros(dftxc.values.shape)

                for d in range(len(dftxc)):
                    # total number of transmissions on index day d
                    n_tx = dftxcv[d]
                    # compute the number that are reported
                    reporting_probability = 0.1
                    n_reported_tx = np.random.binomial(n_tx, reporting_probability)

                    # compute when they are reported
                    delays_days = np.random.lognormal(mean=np.log(8), sigma=np.log(1.35), size=n_reported_tx)
                    assert type(delays_days) is np.ndarray
                    delays_days = np.round(delays_days).astype(int)
                    for delay in delays_days:
                        if d+delay < len(repcases):
                            repcases[d+delay] += 1

                dfrepcases[c] = repcases

                # for d in range(len(dftxc)):
                #     print(dftxcv[d], repcases[d])
                # print(sum(dftxcv), sum(repcases))
                # quit()

            dfrepcases = dfrepcases.set_index('Date')
            dfcmI = dfcmI.set_index('Date')
            dfcmrepcases = dfrepcases.cumsum()
            dfcmrepcases = dfcmrepcases.reset_index()
            fname = join(output_path, 'FL_cumulative_reported_cases_R0_{}_{}_realization_{}.csv'.format(strR0, shortlong, i))
            dfcmrepcases.to_csv(fname, index=False, quoting=csv.QUOTE_NONNUMERIC, date_format="%Y-%m-%d")
        
