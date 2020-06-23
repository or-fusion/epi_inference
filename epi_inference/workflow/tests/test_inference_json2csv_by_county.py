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

class Test_Inference_JSON2CSV_By_County():
    @classmethod
    def setup_class(cls):
        # change to the test directory
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        rundir = os.path.join(thisfiledir, 'inference_json2csv_by_county')
        os.chdir(rundir)

    @classmethod
    def teardown_class(cls):
        # return to the previous directory
        os.chdir(cls._origdir)

    def test_inference_json2csv_by_county(self):
        args = Options()
        args.block = 'all'
        args.config_file = './workflows/inference_json2csv_by_county.yml'
        args.verbose = True
        driver.run(args)

        for f in csv_files():
            output_file = os.path.join('results', f)
            baseline_file = os.path.join('baseline', f) 
            compare_csv(output_file, baseline_file)

def csv_files():
    ret = ['inference_json2csv_by_county/estimated_beta_county_12001.csv',
           'inference_json2csv_by_county/estimated_beta_county_12003.csv',
           'inference_json2csv_by_county/estimated_beta_county_12005.csv',
           'inference_json2csv_by_county/estimated_beta_county_12007.csv',
           'inference_json2csv_by_county/estimated_beta_county_12009.csv',
           'inference_json2csv_by_county/estimated_beta_county_12011.csv',
           'inference_json2csv_by_county/estimated_beta_county_12013.csv',
           'inference_json2csv_by_county/estimated_beta_county_12015.csv',
           'inference_json2csv_by_county/estimated_beta_county_12017.csv',
           'inference_json2csv_by_county/estimated_beta_county_12019.csv',
           'inference_json2csv_by_county/estimated_beta_county_12021.csv',
           'inference_json2csv_by_county/estimated_beta_county_12023.csv',
           'inference_json2csv_by_county/estimated_beta_county_12027.csv',
           'inference_json2csv_by_county/estimated_beta_county_12029.csv',
           'inference_json2csv_by_county/estimated_beta_county_12031.csv',
           'inference_json2csv_by_county/estimated_beta_county_12033.csv',
           'inference_json2csv_by_county/estimated_beta_county_12035.csv',
           'inference_json2csv_by_county/estimated_beta_county_12037.csv',
           'inference_json2csv_by_county/estimated_beta_county_12039.csv',
           'inference_json2csv_by_county/estimated_beta_county_12041.csv',
           'inference_json2csv_by_county/estimated_beta_county_12043.csv',
           'inference_json2csv_by_county/estimated_beta_county_12045.csv',
           'inference_json2csv_by_county/estimated_beta_county_12047.csv',
           'inference_json2csv_by_county/estimated_beta_county_12049.csv',
           'inference_json2csv_by_county/estimated_beta_county_12051.csv',
           'inference_json2csv_by_county/estimated_beta_county_12053.csv',
           'inference_json2csv_by_county/estimated_beta_county_12055.csv',
           'inference_json2csv_by_county/estimated_beta_county_12057.csv',
           'inference_json2csv_by_county/estimated_beta_county_12059.csv',
           'inference_json2csv_by_county/estimated_beta_county_12061.csv',
           'inference_json2csv_by_county/estimated_beta_county_12063.csv',
           'inference_json2csv_by_county/estimated_beta_county_12065.csv',
           'inference_json2csv_by_county/estimated_beta_county_12067.csv',
           'inference_json2csv_by_county/estimated_beta_county_12069.csv',
           'inference_json2csv_by_county/estimated_beta_county_12071.csv',
           'inference_json2csv_by_county/estimated_beta_county_12073.csv',
           'inference_json2csv_by_county/estimated_beta_county_12075.csv',
           'inference_json2csv_by_county/estimated_beta_county_12077.csv',
           'inference_json2csv_by_county/estimated_beta_county_12079.csv',
           'inference_json2csv_by_county/estimated_beta_county_12081.csv',
           'inference_json2csv_by_county/estimated_beta_county_12083.csv',
           'inference_json2csv_by_county/estimated_beta_county_12085.csv',
           'inference_json2csv_by_county/estimated_beta_county_12086.csv',
           'inference_json2csv_by_county/estimated_beta_county_12087.csv',
           'inference_json2csv_by_county/estimated_beta_county_12089.csv',
           'inference_json2csv_by_county/estimated_beta_county_12091.csv',
           'inference_json2csv_by_county/estimated_beta_county_12093.csv',
           'inference_json2csv_by_county/estimated_beta_county_12095.csv',
           'inference_json2csv_by_county/estimated_beta_county_12097.csv',
           'inference_json2csv_by_county/estimated_beta_county_12099.csv',
           'inference_json2csv_by_county/estimated_beta_county_12101.csv',
           'inference_json2csv_by_county/estimated_beta_county_12103.csv',
           'inference_json2csv_by_county/estimated_beta_county_12105.csv',
           'inference_json2csv_by_county/estimated_beta_county_12107.csv',
           'inference_json2csv_by_county/estimated_beta_county_12109.csv',
           'inference_json2csv_by_county/estimated_beta_county_12111.csv',
           'inference_json2csv_by_county/estimated_beta_county_12113.csv',
           'inference_json2csv_by_county/estimated_beta_county_12115.csv',
           'inference_json2csv_by_county/estimated_beta_county_12117.csv',
           'inference_json2csv_by_county/estimated_beta_county_12119.csv',
           'inference_json2csv_by_county/estimated_beta_county_12121.csv',
           'inference_json2csv_by_county/estimated_beta_county_12123.csv',
           'inference_json2csv_by_county/estimated_beta_county_12125.csv',
           'inference_json2csv_by_county/estimated_beta_county_12127.csv',
           'inference_json2csv_by_county/estimated_beta_county_12129.csv',
           'inference_json2csv_by_county/estimated_beta_county_12131.csv',
           'inference_json2csv_by_county/estimated_beta_county_12133.csv']

    return ret

