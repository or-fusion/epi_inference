import pytest
import os
import os.path
import pandas as pd
import shutil
from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options

from epi_inference.engine import driver
from epi_inference.util import compare_csv, compare_json

keepfiles = False

class TestReconstruct():
    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)

    def test_reconstruct_deterministic(self):
        args = Options()
        args.block = 'deterministic'
        args.config_file = './config_files/reconstruct_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the json files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_json('./output/recon_countydata1_all.json', './baseline/recon_countydata1_all.json')
        outputdf, golddf = compare_json('./output/recon_countydata1_12011.json', './baseline/recon_countydata1_12011.json')
    
        # cleanup the files we created
        if not keepfiles:
            os.remove('./output/recon_countydata1_all.json')
            os.remove('./output/recon_countydata1_all_meta.yml')
            os.remove('./output/recon_countydata1_12011.json')
            os.remove('./output/recon_countydata1_12011_meta.yml')

    def test_reconstruct_stochastic(self):
        args = Options()
        args.block = 'stochastic'
        args.config_file = './config_files/reconstruct_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the json files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_json('./output/recon_stoch_countydata1_all_38479387.json', './baseline/recon_stoch_countydata1_all_38479387.json')
        outputdf, golddf = compare_json('./output/recon_stoch_countydata1_all_39847938.json', './baseline/recon_stoch_countydata1_all_39847938.json')
        outputdf, golddf = compare_json('./output/recon_stoch_countydata1_12011.json', './baseline/recon_stoch_countydata1_12011.json')
    
        # cleanup the files we created
        if not keepfiles:
            os.remove('./output/recon_stoch_countydata1_all_38479387.json')
            os.remove('./output/recon_stoch_countydata1_all_38479387_meta.yml')
            os.remove('./output/recon_stoch_countydata1_all_39847938.json')
            os.remove('./output/recon_stoch_countydata1_all_39847938_meta.yml')
            os.remove('./output/recon_stoch_countydata1_12011.json')
            os.remove('./output/recon_stoch_countydata1_12011_meta.yml')

    def test_recon_json2csv(self):
        args = Options()
        args.block = 'json2csv'
        args.config_file = './config_files/reconstruct_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the json files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/recon_stoch_countydata1_all_38479387_flatten.csv', './baseline/recon_stoch_countydata1_all_38479387_flatten.csv')
        outputdf, golddf = compare_csv('./output/recon_stoch_countydata1_all_flatten.csv', './baseline/recon_stoch_countydata1_all_flatten.csv')
        outputdf, golddf = compare_csv('./output/recon_stoch_countydata1_all_38479387_narrow.csv', './baseline/recon_stoch_countydata1_all_38479387_narrow.csv')
        outputdf, golddf = compare_csv('./output/recon_stoch_countydata1_all_narrow.csv', './baseline/recon_stoch_countydata1_all_narrow.csv')
    
        df_flatten = pd.read_csv('./output/recon_stoch_countydata1_all_38479387_flatten.csv')
        df_narrow = pd.read_csv('./output/recon_stoch_countydata1_all_38479387_narrow.csv')

        assert(df_narrow.shape[1] == (df_flatten.shape[1]-8+2))
        assert(df_narrow.shape[0] == 8*df_flatten.shape[0])
        assert(df_narrow.size == (8*df_flatten.shape[0]) * (df_flatten.shape[1]-8+2))

        # cleanup the files we created
        if not keepfiles:
            os.remove('./output/recon_stoch_countydata1_all_38479387_flatten.csv')
            os.remove('./output/recon_stoch_countydata1_all_flatten.csv')
            os.remove('./output/recon_stoch_countydata1_all_38479387_narrow.csv')
            os.remove('./output/recon_stoch_countydata1_all_narrow.csv')
        
    def test_recon_summary(self):
        args = Options()
        args.block = 'recon_summary'
        args.config_file = './config_files/reconstruct_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the json files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/recon_stoch_countydata1_all_summary.csv', './baseline/recon_stoch_countydata1_all_summary.csv')
    
        # cleanup the files we created
        if not keepfiles:
            os.remove('./output/recon_stoch_countydata1_all_summary.csv')

