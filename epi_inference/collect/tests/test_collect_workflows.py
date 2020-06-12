import pytest
import os
import os.path
import pandas as pd
import shutil
from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options
from epi_inference.engine import driver

# collect yaml
# dir, county, days_before_first, days_after_first, output

def compare_csv(output, gold, index_col, check_exact=False):
    outputdf = pd.read_csv(output, index_col=index_col)
    golddf = pd.read_csv(gold, index_col=index_col)
    pd.testing.assert_frame_equal(left=outputdf, right=golddf, check_exact=check_exact)
    return outputdf, golddf


class TestCollect():
    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)
    """
    def setup_method(self, method):
        if os.path.exists('output'):            # pragma: no cover
            shutil.rmtree('output')

    def teardown_method(self, method):
        if os.path.exists('output'):            # pragma: no cover
            shutil.rmtree('output')
    """
    def test_collect_exp(self):
        # run a collection of data for 24031
        args = Options()
        args.block = 'collect'
        args.config_file = './config_files/collect_exp.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/expdata3_1_all.csv', './baseline/expdata3_1_all.csv', index_col='Date')
        assert outputdf.shape[0] == 123
        outputdf, golddf = compare_csv('./output/expdata3_2_all.csv', './baseline/expdata3_2_all.csv', index_col='Date')
        assert outputdf.shape[0] == 123
        outputdf, golddf = compare_csv('./output/expdata3_3_all.csv', './baseline/expdata3_3_all.csv', index_col='Date')
        assert outputdf.shape[0] == 123
        outputdf, golddf = compare_csv('./output/expdata3_1_7before.csv', './baseline/expdata3_1_7before.csv', index_col='Date')
        assert outputdf.shape[0] == 123
        outputdf, golddf = compare_csv('./output/expdata3_1_15after.csv', './baseline/expdata3_1_15after.csv', index_col='Date')
        assert outputdf.shape[0] == 17
        compare_csv('./output/expdata3_1_7before_15after.csv', './baseline/expdata3_1_7before_15after.csv', index_col='Date')
    
        outputdf, golddf = compare_csv('./output/expdata3_1_col3_7before.csv', './baseline/expdata3_1_col3_7before.csv', index_col='Date')
        assert outputdf.shape[0] == 120
        outputdf, golddf = compare_csv('./output/expdata3_1_col3_15after.csv', './baseline/expdata3_1_col3_15after.csv', index_col='Date')
        assert outputdf.shape[0] == 26
        outputdf, golddf = compare_csv('./output/expdata3_1_col3_7before_15after.csv', './baseline/expdata3_1_col3_7before_15after.csv', index_col='Date')
        assert outputdf.shape[0] == 23
    
        # cleanup the files we created
        os.remove('./output/expdata3_1_all.csv')
        os.remove('./output/expdata3_1_all_meta.yml')
        os.remove('./output/expdata3_1_all_geodata.csv')
        os.remove('./output/expdata3_2_all.csv')
        os.remove('./output/expdata3_2_all_meta.yml')
        os.remove('./output/expdata3_2_all_geodata.csv')
        os.remove('./output/expdata3_3_all.csv')
        os.remove('./output/expdata3_3_all_meta.yml')
        os.remove('./output/expdata3_3_all_geodata.csv')
        os.remove('./output/expdata3_1_7before.csv')
        os.remove('./output/expdata3_1_7before_meta.yml')
        os.remove('./output/expdata3_1_7before_geodata.csv')
        os.remove('./output/expdata3_1_15after.csv')
        os.remove('./output/expdata3_1_15after_meta.yml')
        os.remove('./output/expdata3_1_15after_geodata.csv')
        os.remove('./output/expdata3_1_7before_15after.csv')
        os.remove('./output/expdata3_1_7before_15after_meta.yml')
        os.remove('./output/expdata3_1_7before_15after_geodata.csv')
        os.remove('./output/expdata3_1_col3_7before.csv')
        os.remove('./output/expdata3_1_col3_7before_meta.yml')
        os.remove('./output/expdata3_1_col3_7before_geodata.csv')
        os.remove('./output/expdata3_1_col3_15after.csv')
        os.remove('./output/expdata3_1_col3_15after_meta.yml')
        os.remove('./output/expdata3_1_col3_15after_geodata.csv')
        os.remove('./output/expdata3_1_col3_7before_15after.csv')
        os.remove('./output/expdata3_1_col3_7before_15after_meta.yml')
        os.remove('./output/expdata3_1_col3_7before_15after_geodata.csv')
    
        # check that we cleaned everything up
        #files_remaining = [f for f in os.listdir('./output')]
        #assert len(files_remaining) == 0
    
    def test_collect_case(self):
        args = Options()
        args.block = 'collect'
        args.config_file = './config_files/collect_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/countydata1_FL.csv', './baseline/countydata1_FL.csv', index_col='Date')
        assert outputdf.shape[0] == 40

        os.remove('./output/countydata1_FL.csv')
        os.remove('./output/countydata1_FL_meta.yml')

        # check that we cleaned everything up
        #files_remaining = [f for f in os.listdir('./output')]
        #assert len(files_remaining) == 0
