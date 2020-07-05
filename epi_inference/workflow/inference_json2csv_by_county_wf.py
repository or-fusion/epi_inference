import sys
import os.path
try:
    import ujson as json
except:
    import json
import csv
import glob
import pandas as pd
import numpy as np

from epi_inference.engine.task import Task
from epi_inference.engine.task_registry import register_task

def create_inference_csv_by_county(input_json_filespec, output_dir, low_inf_threshold, min_data_days_threshold):
    """
    This function converts json files created by inference of the stochastic
    reconstruction data and produces csv files for each county
    """
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #
    # loop through all the json files and append the necessary fields
    # to a dataframe
    #
    json_filenames = list()
    for filename in glob.glob(input_json_filespec):
        if not os.path.exists(filename):
            raise RuntimeError("ERROR: Inference JSON file does not exist: "+ filename)
        json_filenames.append(filename)
    json_filenames = sorted(json_filenames)

    counties_set = None
    county_dfs = dict()
    for j, jsonfname in enumerate(json_filenames):
        with open(jsonfname, 'r') as fd:
            d = json.load(fd)

        # this is not the best way to do this - we need to have the seed in the json output
        seed = jsonfname[:-5][-7:]
        try:
            test = int(seed)
        except:
            print('inference_json2csv_by_county currently only works with json'
                  ' files produced with inference of stochastic reconstructions'
                  ' since it needs the seed')
            raise
        
        print('... processing seed', seed, j, '/', len(json_filenames))

        # grab a list of counties the first time through
        # to check that each json has the same county list
        if counties_set is None:
            counties_set = set(d.keys())
            sorted_counties = sorted(counties_set)

        # check that all counties are included (and only those counties)
        assert sorted(d.keys()) == sorted_counties

        # build the county dictionaries for this seed
        for c in sorted_counties:
            dc = dict(d[c])
            dc['FIPS'] = [dc['FIPS']]*len(dc['date'])
            dc['window_days'] = [dc['window_days']]*len(dc['date'])
            dc['seed'] = [seed]*len(dc['date'])
            # change the name from beta to raw_est_beta
            dc['raw_est_beta'] = dc['beta']
            del dc['beta']
            # get the dates, changing the name from date to dates
            dc['dates'] = dc['date']
            del dc['date']

            # build the filtered data
            low_inf = list()
            low_data_days = list()
            filtered_est_beta = list()
            for i in range(len(dc['dates'])):
                inf = dc['infections_in_window'][i]
                data_days = dc['days_since_first_reported'][i]
                filtered_beta = dc['raw_est_beta'][i]
                if inf < low_inf_threshold:
                    low_inf.append(True)
                    filtered_beta = None
                else:
                    low_inf.append(False)

                if data_days < min_data_days_threshold:
                    low_data_days.append(True)
                    filtered_beta = None
                else:
                    low_data_days.append(False)

                filtered_est_beta.append(filtered_beta)
            dc['low_infections_in_window(<{})'.format(low_inf_threshold)] = low_inf
            dc['low_data_days(<{}days)'.format(min_data_days_threshold)] = low_data_days
            dc['filtered_est_beta'] = filtered_est_beta

            dfc = pd.DataFrame(dc)
            dfc.fillna(value=np.nan, inplace=True)
            if c not in county_dfs:
                county_dfs[c] = dict()

            county_dfs[c][seed] = dfc


    description_file = os.path.join(output_dir, 'inference_description.txt')
    fd = open(description_file, 'w')
    
    for i,c in enumerate(sorted_counties):
        print('... processing county', c, i, '/', len(sorted_counties))

        # concatenate across all the seeds
        county_df = pd.concat(county_dfs[c].values())

        print('County:', c)
        description_str = county_df.describe(include=None, exclude=[list])
        fd.write(description_str.to_string())
        fd.write('\n')

        fname = os.path.join(output_dir, 'estimated_beta_county_{}.csv'.format(c))
        county_df.to_csv(fname, quoting=csv.QUOTE_NONNUMERIC, index=False)

        # get the status across the seeds
        q05_raw = county_df.groupby('dates')['raw_est_beta'].quantile(.05)
        q05_raw.rename('q05_raw_est_beta', inplace=True)
        q25_raw = county_df.groupby('dates')['raw_est_beta'].quantile(.25)
        q25_raw.rename('q25_raw_est_beta', inplace=True)
        q50_raw = county_df.groupby('dates')['raw_est_beta'].quantile(.50)
        q50_raw.rename('q50_raw_est_beta', inplace=True)
        qmean_raw = county_df.groupby('dates')['raw_est_beta'].mean()
        qmean_raw.rename('qmean_raw_est_beta', inplace=True)
        q75_raw = county_df.groupby('dates')['raw_est_beta'].quantile(.75)
        q75_raw.rename('q75_raw_est_beta', inplace=True)
        q95_raw = county_df.groupby('dates')['raw_est_beta'].quantile(.95)
        q95_raw.rename('q95_raw_est_beta', inplace=True)

        q05_filtered = county_df.groupby('dates')['filtered_est_beta'].quantile(.05)
        q05_filtered.rename('q05_filtered_est_beta', inplace=True)
        q25_filtered = county_df.groupby('dates')['filtered_est_beta'].quantile(.25)
        q25_filtered.rename('q25_filtered_est_beta', inplace=True)
        q50_filtered = county_df.groupby('dates')['filtered_est_beta'].quantile(.50)
        q50_filtered.rename('q50_filtered_est_beta', inplace=True)
        qmean_filtered = county_df.groupby('dates')['filtered_est_beta'].mean()
        qmean_filtered.rename('qmean_filtered_est_beta', inplace=True)
        q75_filtered = county_df.groupby('dates')['filtered_est_beta'].quantile(.75)
        q75_filtered.rename('q75_filtered_est_beta', inplace=True)
        q95_filtered = county_df.groupby('dates')['filtered_est_beta'].quantile(.95)
        q95_filtered.rename('q95_filtered_est_beta', inplace=True)

        county_df_summary = pd.concat([q05_filtered, q25_filtered,
                                       q50_filtered, q75_filtered,
                                       q95_filtered, qmean_filtered,
                                       q05_raw, q25_raw,
                                       q50_raw, q75_raw,
                                       q95_raw, qmean_raw],
                                      axis=1)
        county_df_summary = county_df_summary.reset_index()

        fname = os.path.join(output_dir, 'summary_estimated_beta_county_{}.csv'.format(c))
        county_df_summary.to_csv(fname, quoting=csv.QUOTE_NONNUMERIC, index=False)

    fd.close()


class Inference_JSON2CSV_By_County_Workflow(Task):
    def __init__(self):
        Task.__init__(self, "inference_json2csv_by_county",
            "Create a CSV file of inference results for each county.")

    def validate(self, CONFIG):
        valid_options = set(['input_json', 'output_dir', 'low_infection_threshold', 'factors', 'factor_levels', 'workflow', 'verbose'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        self._warnings = []
        self.validate(CONFIG)
        create_inference_csv_by_county(
            input_json_filespec=CONFIG['input_json'],
            output_dir=CONFIG.get('output_dir', None),
            low_inf_threshold=int(CONFIG.get('low_infection_threshold', 0)),
            min_data_days_threshold=int(CONFIG.get('min_data_days_threshold', 14))
            )

register_task(Inference_JSON2CSV_By_County_Workflow())

if __name__ == '__main__':
    json_filespec = sys.argv[1]
    output_dir = './json2csv-test-results'
    low_infection_threshold = 100
    create_inference_csv_by_county(
        input_json_filespec=json_filespec,
        output_dir=output_dir,
        low_inf_threshold=low_infection_threshold,
        min_data_days_threshold=14
        )
