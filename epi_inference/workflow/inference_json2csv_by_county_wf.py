import sys
import os.path
import json
import csv
import glob
import pandas as pd
import yaml
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
import numpy as np

from epi_inference.engine.task import Task
from epi_inference.engine.task_registry import register_task

def create_inference_csv_by_county(*, input_json_filespec, seeds_yml, output_dir, low_inf_threshold):
    """
    This function converts json files created by inference of 
    reconstruction data and produces csv files for each county
    """
    #
    # Get seeds
    #
    with open(seeds_yml, 'r') as INPUT:
        seeds = yaml.save_load(INPUT)
    #
    # Collect the JSON filenames
    #
    json_filenames = {}
    for seed in seeds:
        filename = input_json_filespec % seed
        if not os.path.exists(filename):
            raise RuntimeError("ERROR: Inference JSON file does not exist: "+ filename)
        json_filenames[seed] = filename
    #
    # loop through all the json files and append the necessary fields
    # to a dataframe
    #
    counties_set = None
    county_dfs = dict()
    j = 0
    for seed in sorted(seeds):
        j += 1
        with open(json_filenames[seed], 'r') as fd:
            d = json.load(fd)

        # this is not the best way to do this - we need to have the seed in the json output
        seed = int(seed)
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

            # build the low infections filter
            low_inf = list()
            filtered_est_beta = list()
            
            for i in range(len(dc['dates'])):
                inf = dc['infections_in_window'][i]
                beta = dc['raw_est_beta'][i]
                if inf < low_inf_threshold:
                    low_inf.append(True)
                    filtered_est_beta.append(None)
                else:
                    low_inf.append(False)
                    filtered_est_beta.append(beta)
            dc['low_infections_in_window(<{})'.format(low_inf_threshold)] = low_inf
            dc['filtered_est_beta'] = filtered_est_beta

            dfc = pd.DataFrame(dc)
            dfc.fillna(value=np.nan, inplace=True)
            if c not in county_dfs:
                county_dfs[c] = dict()
                
            county_dfs[c][seed] = dfc


    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    summary_file = os.path.join(output_dir, 'inference_summary.txt')
    fd = open(summary_file, 'w')
    
    for i,c in enumerate(sorted_counties):
        print('... processing county', c, i, '/', len(sorted_counties))

        # concatenate across all the seeds
        county_df = pd.concat(county_dfs[c].values())

        print('County:', c)
        summary_str = county_df.describe(include='all')
        fd.write(summary_str.to_string())
        fd.write('\n')

        fname = os.path.join(output_dir, 'estimated_beta_county_{}.csv'.format(c))
        county_df.to_csv(fname, quoting=csv.QUOTE_NONNUMERIC, index=False)
    fd.close()


class Inference_JSON2CSV_By_County_Workflow(Task):
    def __init__(self):
        Task.__init__(self, "inference_json2csv_by_county",
            "Create a CSV file of inference results for each county.")

    def validate(self, CONFIG):
        valid_options = set(['input_json', 'seeds_yml', 'output_dir', 'low_infection_threshold', 'factors', 'factor_levels', 'workflow', 'verbose'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        create_inference_csv_by_county(
            input_json_filespec=CONFIG['input_json'],
            seeds_yml=CONFIG['seeds_yml'],
            output_dir=CONFIG.get('output_dir', None),
            low_inf_threshold=int(CONFIG.get('low_infection_threshold', 0))
            )

register_task(Inference_JSON2CSV_By_County_Workflow())

if __name__ == '__main__':
    json_filespec = sys.argv[1]
    output_dir = './json2csv-test-results'
    low_infection_threshold = 100
    create_inference_csv_by_county(
        input_json_filespec=json_filespec,
        output_dir=output_dir,
        low_inf_threshold=low_infection_threshold
        )
