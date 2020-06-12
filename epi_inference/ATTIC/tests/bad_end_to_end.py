import pytest
import subprocess
import os
import os.path
import json
import pandas as pd
import numpy as np

from pyomo.common import fileutils as fileutils

def _cleanup_output(files):
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

class TestEndToEnd():

    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)

    def test_inference_simulated_county_different_beta(self):
        # list of expected output files for cleanup
        expected_output_files = ['./output/results_inference_simulated_county_different_beta.csv',
                        './output/results_inference_simulated_county_different_beta_meta.yml']

        # cleanup any output that may be lingering from old runs
        _cleanup_output(expected_output_files)

        # execute epiinf
        cmd = ['epiinf', 'inference', './config_files/inference_simulated_county_different_beta.yml']
        subprocess.run(cmd)

        # read the output that was produced
        assert os.path.isfile('./output/results_inference_simulated_county_different_beta.csv')
        df = pd.read_csv('./output/results_inference_simulated_county_different_beta.csv', dtype={"FIPS":'str'})

        # expected output
        expected_df = {"est_beta": [0.25, 0.5, 0.75, 1.0, 1.25, 1.5],
                       "status": ['ok' for i in range(6)],
                       "population": [1000000]*6,
                       "total_cases": [2.579078, 6.097038, 12.661745, 23.655522, 40.762154, 65.994705],
                       "FIPS": ['025', '050', '075', '100', '125', '150']}
        expected_df = pd.DataFrame.from_dict(expected_df)

        # compare the estimated beta values
        pd.testing.assert_series_equal(left=df['est_beta'], right=expected_df['est_beta'], check_exact=False)
        pd.testing.assert_series_equal(left=df['total_cases'], right=expected_df['total_cases'], check_exact=False)
        pd.testing.assert_series_equal(left=df['status'], right=expected_df['status'], check_exact=False)

        # cleanup the output 
        _cleanup_output(expected_output_files)

    def test_inference_simulated_county_different_beta_int(self):
        # list of expected output files for cleanup
        expected_output_files = ['./output/results_inference_simulated_county_different_beta_int.csv',
                        './output/results_inference_simulated_county_different_beta_int_meta.yml']

        # cleanup any output that may be lingering from old runs
        _cleanup_output(expected_output_files)

        # execute epiinf
        cmd = ['epiinf', 'inference', './config_files/inference_simulated_county_different_beta_int.yml']
        subprocess.run(cmd)

        # read the output that was produced
        assert os.path.isfile('./output/results_inference_simulated_county_different_beta_int.csv')
        df = pd.read_csv('./output/results_inference_simulated_county_different_beta_int.csv', dtype={"FIPS":'str'})

        # expected output
        expected_df = {"est_beta": [0.256255, 0.507560, 0.802103, 1.015384, 1.271785, 1.480856],
                       "status": ['ok' for i in range(6)],
                       "population": [1000000]*6,
                       "total_cases": [3.0, 6.0, 13.0, 24.0, 41.0, 66.0],
                       "FIPS": ['025', '050', '075', '100', '125', '150']}
        expected_df = pd.DataFrame.from_dict(expected_df)

        # compare the produced output with the expected output
        pd.testing.assert_series_equal(left=df['est_beta'], right=expected_df['est_beta'], check_exact=False)
        pd.testing.assert_series_equal(left=df['total_cases'], right=expected_df['total_cases'], check_exact=False)
        pd.testing.assert_series_equal(left=df['status'], right=expected_df['status'], check_exact=False)

        # cleanup the output 
        _cleanup_output(expected_output_files)

    def test_inference_florida_mobility(self):
        # list of expected output files for cleanup
        expected_output_files = ['./output/florida_inference_mobility_2.json',
                        './output/florida_inference_mobility_2_meta.yml']

        # cleanup any output that may be lingering from old runs
        _cleanup_output(expected_output_files)

        # execute epiinf
        cmd = ['epiinf', 'inference', './config_files/florida_inference_mobility.yml']
        subprocess.run(cmd)

        # check that the output files were produced
        for f in expected_output_files:
            assert os.path.isfile(f)

        # read the json output
        with open('./output/florida_inference_mobility_2.json', 'r') as fd:
            data = json.load(fd)

        # loop over each run (different R0)
        expected_avg_beta = [1.36/4, 1.51/4, 1.81/4, 2.02/4, 2.17/4, 2.37/4, 2.55/4, 2.70/4]
        betas_by_r0 = [list()]*len(expected_avg_beta)
        for r, rd in enumerate(data):
            r0 = float(rd['R0'])
            for d in rd['results']:
                if d['total_cases'] > 20:
                    betas_by_r0[r].append(d['est_beta'])
            avg_beta = sum(betas_by_r0[r])/len(betas_by_r0[r])
            print(avg_beta*4, r0)
            assert abs(expected_avg_beta[r] - avg_beta) <= 0.01*expected_avg_beta[r]
            assert abs(float(r0)-avg_beta*4) <= 0.1*r0
