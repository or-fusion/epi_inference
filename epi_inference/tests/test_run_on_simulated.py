import pytest
import os
import os.path

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
        # TODO: Add automated comparison of results with simulation parameters
        # If these tests pass, it is safe to update the gold standards below
        
        #
        # compare all the json files
        #
        for f in results_json_files():
            output_file = os.path.join('results', f)
            baseline_file = os.path.join('gold_results', f) 
            compare_json(output_file, baseline_file, abs_tol=1e-6)

        #
        # compare all the csv files
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
