import pytest
import os
import os.path
try:
    import ujson as json
except:
    import json
import yaml
import numpy as np
import pandas as pd

from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options

from epi_inference.engine import driver
from epi_inference.util import compare_json, compare_csv

skip_new_files = False

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

        # remove old results files used in comparison
        res_files = _walk_files('results/run_on_simulated', '.json')
        for f in res_files:
            os.remove(os.path.join('results/run_on_simulated', f))
        res_files = _walk_files('results/run_on_simulated', '.csv')
        for f in res_files:
            os.remove(os.path.join('results/run_on_simulated', f))

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
            with open('./results/run_on_simulated/reconstruct_stochastic_confirmed/recon_stochastic_10_{}.json'.format(str(sd)),'r') as fd:
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
                            # print('Outside confidence:', cname, sname, dt, v, ':', ql, mn, qu)
                        if abs(v-mn) > 10 and abs(v-mn)/max(1,v) > 0.2:
                            mean_error_count += 1
                            # print('Error with mean > 20%:', cname, sname, dt, v, ':', ql, mn, qu)

        print('Reconstruction fraction outside 95%:', outside_count/total_count)
        print('Reconstruction fraction with error in mean > 20%:', mean_error_count/total_count)
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
            with open('./results/run_on_simulated/inference_stochastic_confirmed/inference_mobility_window_10_{}.json'.format(str(sd)),'r') as fd:
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
                        # print('Beta outside confidence:', cname, 'beta', dt, v, ':', ql, mn, qu, '# Not None:', len(data) )
                    if abs(v-mn)/max(1,v) > 0.2:
                        mean_error_count += 1
                        # print('Beta error with mean > 20%:', cname, 'beta', dt, v, ':', ql, mn, qu, '# Not None:', len(data))

        print('Inference fraction outside 95%:', outside_count/total_count)
        print('Inference fraction with error in mean > 20%:', mean_error_count/total_count)
        assert outside_count / total_count < 0.05
        assert mean_error_count / total_count < 0.15

        #
        # compare all the json files (gold standard)
        #
        res_json_files = _walk_files('results/run_on_simulated', '.json')
        baseline_json_files = set(_walk_files('baseline/run_on_simulated', '.json'))
        for f in res_json_files:
            if not skip_new_files:
                assert f in baseline_json_files # if this fails then there are files in the new results that are not in the baseline

            if f in baseline_json_files:
                baseline_json_files.remove(f)
                res_file = os.path.join('results/run_on_simulated', f)
                baseline_file = os.path.join('baseline/run_on_simulated', f)
                compare_json(res_file, baseline_file, abs_tol=1e-6)

        if len(baseline_json_files) != 0:
            print(baseline_json_files)
        assert len(baseline_json_files) == 0 # if this fails, then there are files in the baseline that did not appear in the new results

        #
        # compare all the csv files (gold standard)
        #
        res_csv_files = _walk_files('results/run_on_simulated', '.csv')
        baseline_csv_files = set(_walk_files('baseline/run_on_simulated', '.csv'))
        for f in res_csv_files:
            if not skip_new_files:
                assert f in baseline_csv_files # if this fails then there are files in the new results that are not in the baseline
            if f in baseline_csv_files:
                baseline_csv_files.remove(f)
                res_file = os.path.join('results/run_on_simulated', f)
                baseline_file = os.path.join('baseline/run_on_simulated', f)
                compare_csv(res_file, baseline_file, check_exact=False)

        if len(baseline_csv_files) != 0:
            print(baseline_csv_files)
        assert len(baseline_csv_files) == 0 # if this fails, then there are files in the baseline that did not appear in the new results


    @pytest.mark.skip('skipped - not testing anything that run_on_simulated does not')
    def test_run_on_simulated_half(self):
        args = Options()
        args.block = 'all'
        args.config_file = './workflows/run_on_simulated_half.yml'
        args.verbose = True

        # remove old results files used in comparison
        res_files = _walk_files('results/run_on_simulated_half', '.json')
        for f in res_files:
            os.remove(os.path.join('results/run_on_simulated_half', f))
        res_files = _walk_files('results/run_on_simulated_half', '.csv')
        for f in res_files:
            os.remove(os.path.join('results/run_on_simulated_half', f))

        driver.run(args)

        #
        # compare all the json files (gold standard)
        #
        res_json_files = _walk_files('results/run_on_simulated_half', '.json')
        baseline_json_files = set(_walk_files('baseline/run_on_simulated_half', '.json'))
        for f in res_json_files:
            if not skip_new_files:
                assert f in baseline_json_files # if this fails then there are files in the new results that are not in the baseline

            if f in baseline_json_files:
                baseline_json_files.remove(f)
                res_file = os.path.join('results/run_on_simulated_half', f)
                baseline_file = os.path.join('baseline/run_on_simulated_half', f)
                compare_json(res_file, baseline_file, abs_tol=1e-6)

        assert len(baseline_json_files) == 0 # if this fails, then there are files in the baseline that did not appear in the new results

        #
        # compare all the csv files (gold standard)
        #
        res_csv_files = _walk_files('results/run_on_simulated_half', '.csv')
        baseline_csv_files = set(_walk_files('baseline/run_on_simulated_half', '.csv'))
        for f in res_csv_files:
            if not skip_new_files:
                assert f in baseline_csv_files # if this fails then there are files in the new results that are not in the baseline
            if f in baseline_csv_files:
                baseline_csv_files.remove(f)
                res_file = os.path.join('results/run_on_simulated_half', f)
                baseline_file = os.path.join('baseline/run_on_simulated_half', f)
                compare_csv(res_file, baseline_file, check_exact=False)
                
        assert len(baseline_csv_files) == 0 # if this fails, then there are files in the baseline that did not appear in the new results

def _walk_files(basepath, extension):
    # get all the files of a particular extension in the directory tree
    ret_files = list()

    for path, folders, files in os.walk(basepath):
        for f in files:
            fname, fext = os.path.splitext(f)
            if fext == extension:
                relfname = os.path.relpath(os.path.join(path,f), basepath)
                ret_files.append(relfname)
    return ret_files
