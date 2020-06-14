import pytest
import os
import os.path
import shutil

from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options
from epi_inference.engine import driver
from epi_inference.util import compare_csv, compare_json



class TestInference():
    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)

    def test_mobility_window(self):
        # run a collection of data for 24031
        args = Options()
        args.block = 'mobility_windows'
        args.config_file = './config_files/tests1.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        compare_json('./output/tests1_inference_unsampled_countydata1_all.json', './baseline/tests1_inference_unsampled_countydata1_all.json')
    
        # cleanup the files we created
        os.remove('./output/tests1_inference_unsampled_countydata1_all.json')
        os.remove('./output/tests1_inference_unsampled_countydata1_all_meta.yml')

    def test_select_mobility_window(self):
        # run a collection of data for 24031
        args = Options()
        args.block = 'select_mobility_windows'
        args.config_file = './config_files/tests1.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        compare_json('./output/tests1_inference_unsampled_countydata1_all_select_last.json', './baseline/tests1_inference_unsampled_countydata1_all_select_last.json')
        compare_json('./output/tests1_inference_unsampled_countydata1_all_select_3.json', './baseline/tests1_inference_unsampled_countydata1_all_select_last.json')
    
        # cleanup the files we created
        os.remove('./output/tests1_inference_unsampled_countydata1_all_select_last.json')
        os.remove('./output/tests1_inference_unsampled_countydata1_all_select_last_meta.yml')
        os.remove('./output/tests1_inference_unsampled_countydata1_all_select_3.json')
        os.remove('./output/tests1_inference_unsampled_countydata1_all_select_3_meta.yml')


