import pytest
import os
import os.path
import json

from jsondiff import diff
from pyomo.common import fileutils as fileutils
from pyutilib.misc import Options as Options
from epi_inference.engine import driver
from epi_inference.engine.task import Task
from epi_inference.engine.task_registry import register_task


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


class Task1(Task):

    def __init__(self):
        Task.__init__(self, 'task1', "A task used for testing")

    def run(self, data, args):
        with open(args['output'], "w") as OUTPUT:
            json.dump(args, OUTPUT)

register_task(Task1())


class TestEngine():
    @classmethod
    def setup_class(cls):
        cls._origdir = os.getcwd()
        thisfiledir = fileutils.this_file_dir()
        os.chdir(thisfiledir)

    @classmethod
    def teardown_class(cls):
        os.chdir(cls._origdir)

    def test_test1(self):
        args = Options()
        args.block = 'test1'
        args.config_file = 'tests1.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test1_output1.json', 'test1_baseline1.json')
        os.remove('test1_output1.json')

        compare_json('test1_output2.json', 'test1_baseline2.json')
        os.remove('test1_output2.json')

    def test_test2(self):
        args = Options()
        args.block = 'test2'
        args.config_file = 'tests1.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test2_output1.json', 'test2_baseline1.json')
        os.remove('test2_output1.json')

    def test_test3(self):
        args = Options()
        args.block = 'test3'
        args.config_file = 'tests1.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test1_output1.json', 'test1_baseline1.json')
        os.remove('test1_output1.json')

        compare_json('test1_output2.json', 'test1_baseline2.json')
        os.remove('test1_output2.json')

        compare_json('test2_output1.json', 'test2_baseline1.json')
        os.remove('test2_output1.json')

    def test_test4(self):
        args = Options()
        args.block = 'test4'
        args.config_file = 'tests1.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test1_output1.json', 'test1_baseline1.json')
        os.remove('test1_output1.json')

        compare_json('test1_output2.json', 'test1_baseline2.json')
        os.remove('test1_output2.json')

        compare_json('test2_output1.json', 'test2_baseline1.json')
        os.remove('test2_output1.json')

    def test_test5(self):
        args = Options()
        args.block = 'test5'
        args.config_file = 'tests1.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test1_output1.json', 'test1_baseline1.json')
        os.remove('test1_output1.json')

        compare_json('test1_output2.json', 'test1_baseline2.json')
        os.remove('test1_output2.json')

        compare_json('test2_output1.json', 'test2_baseline1.json')
        os.remove('test2_output1.json')

    def test_test6(self):
        args = Options()
        args.block = 'test6'
        args.config_file = 'tests2.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test6_output1.json', 'test6_baseline1.json')
        os.remove('test6_output1.json')

        compare_json('test6_output2.json', 'test6_baseline2.json')
        os.remove('test6_output2.json')

    def test_test7(self):
        args = Options()
        args.block = 'test7'
        args.config_file = 'tests3.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test7_output1.json', 'test7_baseline1.json')
        os.remove('test7_output1.json')

    def test_test8(self):
        args = Options()
        args.block = 'test8'
        args.config_file = 'tests4.yml'
        args.verbose = True
        driver.run(args)

        compare_json('test8_output1.json', 'test8_baseline1.json')
        os.remove('test8_output1.json')

