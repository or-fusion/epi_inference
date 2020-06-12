import pytest
import os
import os.path
import pandas as pd
import shutil
from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options
from epi_inference.engine import driver
import json
from jsondiff import diff

def compare_csv(output, gold, index_col=None, check_exact=False, sort=True):
    if index_col is None:
        outputdf = pd.read_csv(output)
        golddf = pd.read_csv(gold)
    else:
        outputdf = pd.read_csv(output, index_col=index_col)
        golddf = pd.read_csv(gold, index_col=index_col)

    # the dataframes may be the same, but just in a different order
    if sort:
        columns = list(outputdf.columns)
        outputdf.sort_values(by=columns, inplace=True, ignore_index=True)
        golddf.sort_values(by=columns, inplace=True, ignore_index=True)
    pd.testing.assert_frame_equal(left=outputdf, right=golddf, check_exact=check_exact)
    return outputdf, golddf

def compare_json(output, gold, check_exact=False):            # pragma: no cover
    with open(output,'r') as INPUT:
        outputdf = json.load(INPUT)
    with open(gold,'r') as INPUT:
        golddf = json.load(INPUT)
    d = diff(outputdf, golddf)
    if len(d) != 0:
        print('DIFFERENCES IN JSON')
        print(d)
    assert(len(d) == 0)
    return outputdf, golddf


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

        assert(df_narrow.shape[1] == (df_flatten.shape[1]-7+2))
        assert(df_narrow.shape[0] == 7*df_flatten.shape[0])
        assert(df_narrow.size == (7*df_flatten.shape[0]) * (df_flatten.shape[1]-7+2))

        # cleanup the files we created
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
        os.remove('./output/recon_stoch_countydata1_all_summary.csv')

