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

def compare_csv(output, gold, index_col, cols_to_compare=None, check_exact=False):            # pragma: no cover
    if index_col is None:
        outputdf = pd.read_csv(output)
        golddf = pd.read_csv(gold)
    else:
        outputdf = pd.read_csv(output, index_col=index_col)
        golddf = pd.read_csv(gold, index_col=index_col)

    #print(outputdf)
    #print(golddf)

    if cols_to_compare is None:
        pd.testing.assert_frame_equal(left=outputdf, right=golddf, check_exact=check_exact)
    else:
        for c in cols_to_compare:
            pd.testing.assert_series_equal(left=outputdf[c], right=golddf[c], check_exact=check_exact)

    return outputdf, golddf


class TestInference():
    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)


    @pytest.mark.skip('inconsistencies in simulated data seems to be causing more cases than population')
    def test_inference_decay_lsq_filter(self):
        args = Options()
        args.config_file = './config_files/inference_decay_lsq_filter.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/results_expdata3_1_all.csv', './baseline/results_expdata3_1_all.csv', index_col=None)
        assert outputdf.shape[0] == 10
        outputdf, golddf = compare_csv('./output/results_expdata3_1_all_lastdays.csv', './baseline/results_expdata3_1_all_lastdays.csv', index_col=None)
        assert outputdf.shape[0] == 10
    
        # cleanup the files we created
        os.remove('./output/results_expdata3_1_all.csv')
        os.remove('./output/results_expdata3_1_all_meta.yml')
        os.remove('./output/results_expdata3_1_all_lastdays.csv')
        os.remove('./output/results_expdata3_1_all_lastdays_meta.yml')


    @pytest.mark.skip('inconsistencies in simulated data seems to be causing more cases than population')
    def test_inference_exp(self):
        args = Options()
        args.config_file = './config_files/inference_exp.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/results_expdata3_decay-lsq_all.csv', './baseline/results_expdata3_decay-lsq_all.csv', index_col=None)
        assert outputdf.shape[0] == 10
        outputdf, golddf = compare_csv('./output/results_expdata3_inference_all.csv', './baseline/results_expdata3_inference_all.csv', index_col=None)
        assert outputdf.shape[0] == 20
    
        # cleanup the files we created
        os.remove('./output/results_expdata3_decay-lsq_all.csv')
        os.remove('./output/results_expdata3_decay-lsq_all_meta.yml')
        os.remove('./output/results_expdata3_inference_all.csv')
        os.remove('./output/results_expdata3_inference_all_meta.yml')
    

    def test_inference_case(self):
        args = Options()
        args.config_file = './config_files/inference_case.yml'
        args.verbose = True
        driver.run(args)
    
        # check that the csv files load into dataframes that have the correct numbers and shapes
        outputdf, golddf = compare_csv('./output/results_countydata1_12121_decay-lsq_all.csv', './baseline/results_countydata1_12121_decay-lsq_all.csv',
                                       cols_to_compare=['est_beta', 'status', 'FIPS'], index_col=None)
        assert outputdf.shape[0] == 7
        outputdf, golddf = compare_csv('./output/results_county1_inference_all.csv', './baseline/results_county1_inference_all.csv',
                                       cols_to_compare=['est_beta', 'status', 'FIPS'], index_col=None)
        assert outputdf.shape[0] == 14
        outputdf, golddf = compare_csv('./output/results_countydata1_multinode_all.csv', './baseline/results_countydata1_multinode_all.csv',
                                       cols_to_compare=['est_beta', 'status'], index_col=None)
        assert outputdf.shape[0] == 1
    
        # cleanup the files we created
        os.remove('./output/results_countydata1_12121_decay-lsq_all.csv')
        os.remove('./output/results_countydata1_12121_decay-lsq_all_meta.yml')
        os.remove('./output/results_county1_inference_all.csv')
        os.remove('./output/results_county1_inference_all_meta.yml')
        os.remove('./output/results_countydata1_multinode_all.csv')
        os.remove('./output/results_countydata1_multinode_all_meta.yml')

