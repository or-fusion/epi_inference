import pytest
import os
import os.path
import json
import yaml
import numpy as np
import pandas as pd

from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options

from epi_inference.engine import driver
from epi_inference.util import compare_json, compare_csv

class TestRunOnSimulated():
    @classmethod
    def setup_class(cls):
        # change to the test directory
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        rundir = os.path.join(thisfiledir, 'run_on_simulated')
        os.chdir(rundir)

    @classmethod
    def teardown_class(cls):
        # return to the previous directory
        os.chdir(cls._origdir)

    def test_run_on_simulated(self):
        args = Options()
        args.block = 'all'
        args.config_file = './workflows/run_on_simulated.yml'
        args.verbose = True
        driver.run(args)

        #
        # compare the reconstruction results
        #
        # read in the simulation data
        simdf = pd.read_csv('./simulated_data/SEIIIR_R0_2.25_short_realization_4.csv')
        
        # read in the seeds to help with comparisons
        with open('./config/seeds_50.yml', 'r') as fd:
            seeds = yaml.safe_load(fd)
        seeds = [str(sd) for sd in seeds]

        #
        # Comparison of results with simulation parameters
        # If these tests pass, it is safe to update the gold standards compared later
        series = ['transmissions', 'S', 'E', 'I1', 'I2', 'I3', 'R']

        # read all the json files and compute statistics across the seeds
        recon = dict()
        for sd in seeds:
            with open('./results/run_on_simulated/recon/recon_stochastic_10_{}.json'.format(str(sd)),'r') as fd:
                recon[sd] = json.load(fd)

        counties = list(recon[seeds[0]].keys())
        dates = list(recon[seeds[0]][counties[0]]["dates"])

        # check that the counties are the same across all realizations
        for sd in seeds:
            assert counties == list(recon[sd].keys())

        # check that the dates are the same across all counties, realizations
        for cname in counties:
            for sd in seeds:
                assert dates == recon[sd][cname]['dates']

        # check each county, date, series
        total_count = 0
        outside_count = 0
        mean_error_count = 0
        for cname in counties:
            for sname in series:
                sim = simdf[simdf['comp']==sname]
                for i,dt in enumerate(dates):
                    if dt in sim.Date.values:
                        simdt = sim[sim['Date']==dt]
                        data = list()
                        for sd in seeds:
                            data.append(recon[sd][cname][sname][i])
                        data = np.asarray(data)
                        mn = np.mean(data)
                        ql = np.quantile(data,q=0.025)
                        qu = np.quantile(data,q=0.975)
                        v = float(simdt[cname].values)
                        total_count += 1
                        if (v < ql or v > qu):
                            outside_count += 1
                            print('Outside confidence:', cname, sname, dt, v, ':', ql, mn, qu)
                        if abs(v-mn) > 10 and abs(v-mn)/max(1,v) > 0.2:
                            mean_error_count += 1
                            print('Error with mean > 20%:', cname, sname, dt, v, ':', ql, mn, qu)

        print('Fraction outside 95%:', outside_count/total_count)
        print('Fraction with error in mean > 20%:', mean_error_count/total_count)
        assert outside_count / total_count < 0.05
        assert mean_error_count / total_count < 0.05

        #
        # compare the inference results
        # This comparison is filtered - window must have at least 50 cases,
        # and there must be at least 25 (half) realizations with estimates
        #
        # read all the json files and compute statistics across the seeds
        inference = dict()
        for sd in seeds:
            with open('./results/run_on_simulated/inference/inference_mobility_window_10_{}.json'.format(str(sd)),'r') as fd:
                inference[sd] = json.load(fd)

        counties = list(inference[seeds[0]].keys())
        dates = list(inference[seeds[0]][counties[0]]["date"])

        # check that the counties are the same across all realizations
        for sd in seeds:
            assert counties == list(inference[sd].keys())

        # check that the dates are the same across all counties, realizations
        for cname in counties:
            for sd in seeds:
                assert dates == inference[sd][cname]['date']

        # check each county, date
        total_count = 0
        outside_count = 0
        mean_error_count = 0
        for cname in counties:
            for i,dt in enumerate(dates):
                data = list()
                for sd in seeds:
                    estbeta = inference[sd][cname]['beta'][i]
                    if estbeta is not None and inference[sd][cname]['infections_in_window'][i] > 50: # filter out insufficient data - 100 means there were only ~10 reported cases in window
                        data.append(estbeta)
                if len(data) >= 25: ### # need at least 10 that are not None to get some statistics
                    data = np.asarray(data)
                    mn = np.mean(data)
                    ql = np.quantile(data,q=0.025)
                    qu = np.quantile(data,q=0.975)
                    v = 2.25/4.3
                    total_count += 1
                    if (v < ql or v > qu):
                        outside_count += 1
                        print('Beta outside confidence:', cname, 'beta', dt, v, ':', ql, mn, qu, '# Not None:', len(data) )
                    if abs(v-mn)/max(1,v) > 0.2:
                        mean_error_count += 1
                        print('Beta error with mean > 20%:', cname, 'beta', dt, v, ':', ql, mn, qu, '# Not None:', len(data))

        print('Inference fraction outside 95%:', outside_count/total_count)
        print('Inference fraction with error in mean > 20%:', mean_error_count/total_count)
        assert outside_count / total_count < 0.05
        assert mean_error_count / total_count < 0.15

        #
        # compare all the json files (gold standard)
        #
        for f in results_json_files():
            output_file = os.path.join('results', f)
            baseline_file = os.path.join('gold_results', f) 
            compare_json(output_file, baseline_file, abs_tol=1e-6)

        #
        # compare all the csv files (gold standard)
        #
        for f in results_csv_files():
            output_file = os.path.join('results', f)
            baseline_file = os.path.join('gold_results', f) 
            compare_csv(output_file, baseline_file)

def results_json_files():
    return ['run_on_simulated/recon/recon_deterministic_delay_10.json',
            'run_on_simulated/recon/recon_stochastic_10_1236482.json',
            'run_on_simulated/recon/recon_stochastic_10_2800890.json',
            'run_on_simulated/recon/recon_stochastic_10_2840373.json',
            'run_on_simulated/recon/recon_stochastic_10_3116341.json',
            'run_on_simulated/recon/recon_stochastic_10_4995435.json',
            'run_on_simulated/recon/recon_stochastic_10_6742372.json',
            'run_on_simulated/recon/recon_stochastic_10_7549469.json',
            'run_on_simulated/recon/recon_stochastic_10_8334300.json',
            'run_on_simulated/recon/recon_stochastic_10_8573365.json',
            'run_on_simulated/recon/recon_stochastic_10_9797259.json',
            'run_on_simulated/inference/inference_mobility_window_10_1236482.json',
            'run_on_simulated/inference/inference_mobility_window_10_2800890.json',
            'run_on_simulated/inference/inference_mobility_window_10_2840373.json',
            'run_on_simulated/inference/inference_mobility_window_10_3116341.json',
            'run_on_simulated/inference/inference_mobility_window_10_4995435.json',
            'run_on_simulated/inference/inference_mobility_window_10_6742372.json',
            'run_on_simulated/inference/inference_mobility_window_10_7549469.json',
            'run_on_simulated/inference/inference_mobility_window_10_8334300.json',
            'run_on_simulated/inference/inference_mobility_window_10_8573365.json',
            'run_on_simulated/inference/inference_mobility_window_10_9797259.json',
            'run_on_simulated/inference/inference_mobility_window_unsampled_10.json']

def results_csv_files():
    return ['run_on_simulated/recon/recon_stochastic_summary_10.csv',
            'run_on_simulated/recon/recon_deterministic_delay_10.csv',
            'run_on_simulated/recon/recon_stochastic_10.csv',
            'run_on_simulated/inference/inference_mobility_window_10_1236482.csv',
            'run_on_simulated/inference/inference_mobility_window_10_2800890.csv',
            'run_on_simulated/inference/inference_mobility_window_10_2840373.csv',
            'run_on_simulated/inference/inference_mobility_window_10_3116341.csv',
            'run_on_simulated/inference/inference_mobility_window_10_4995435.csv',
            'run_on_simulated/inference/inference_mobility_window_10_6742372.csv',
            'run_on_simulated/inference/inference_mobility_window_10_7549469.csv',
            'run_on_simulated/inference/inference_mobility_window_10_8334300.csv',
            'run_on_simulated/inference/inference_mobility_window_10_8573365.csv',
            'run_on_simulated/inference/inference_mobility_window_10_9797259.csv']
