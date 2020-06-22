__all__ = ['driver', 'run', 'run_block_multiworkflow']

import os
import sys
import datetime
import argparse
import yaml
import pprint
from pyutilib.misc import Options, timing
import pyutilib.services
try:
    import joblib
    joblib_available = True
except:
    joblib_available = False

from .task_registry import registered_tasks
from .util import factorial_iterator, load_configfile
from .config_parameters import set_config_parameters, get_config_parameters
from .misc import save_metadata

global_tempdir = None

def driver(tmpdir=None, command='unknown'):
    if os.environ.get('TMPDIR',tmpdir) is not None:
        _tmpdir = os.environ.get('TMPDIR',tmpdir)
        if not os.path.exists(_tmpdir):
            os.makedirs(_tmpdir)
        #os.environ['TMPDIR'] = tmpdir
        pyutilib.services.TempfileManager.tempdir = _tmpdir
        global global_tempdir
        global_tempdir = _tmpdir

    Parser = argparse.ArgumentParser(description='inference models')
    Parser.add_argument('-c', '--catch-errors', action='store_true', default=False,
                    help='Catch exceptions')
    Parser.add_argument('-q', '--quiet', action='store_true', default=False,
                    help='Suppress diagnostic')
    Parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Verbosity flag')
    Parser.add_argument('-P', '--parallelize-workflows', action='store_true', default=False,
                    help='Parallelize the execution of all workflows.  This option parallelizes all workflows and their factors.')
    Parser.add_argument('-p', '--parallelize-factors', action='store_true', default=False,
                    help='Parallelize the execution of factors for each workflow if there are enough of them (>= MIN_PARALLEL_WORKFLOWS).')
    Parser.add_argument('--np', action='store', type=int, default=2,
                    help='Number of processors used to parallelize workflows. (default=2)')
    Parser.add_argument('--min-parallel-workflows', action='store', type=int, default=2,
                    help='Minimum number of workflows that are needed to use parallelization. (default=2)')
    Parser.add_argument('--help-workflows', action='store_true', default=False,
                    help='Print all available workflows')
    Parser.add_argument('block', help='Name of block in the config file that will be executed', default=None)
    Parser.add_argument('config_file', help='YAML configuration file', default=None)
    #
    # Parse sys.argv
    #
    args = Parser.parse_args()
    if not args.quiet:
        timing.tic('Starting '+command)
    if args.catch_errors:
        run(args)
    else:
        try:
            run(args)
        except Exception as e:
            print("Error: "+str(e))
    if not args.quiet:
        timing.toc('Finishing '+command)


def run_workflow(task, CONFIG, data, tempdir):
    if tempdir:
        pyutilib.services.TempfileManager.tempdir = tempdir
    print("\nExecuting Workflow: "+task.name)
    task.run( data, CONFIG )


def run_block_serial(workflow_blocks, args, config_parameters):
    #
    # Serially execute all workflows and their factors
    #
    tasks = registered_tasks()
    for cargs in workflow_blocks[args.block]:
        workflow = cargs.get('workflow', None)
        if workflow is None:
            print("WARNING: No workflow specified in '%s' block" % args.block)
            continue
        if not 'verbose' in cargs:
            cargs['verbose'] = args.verbose
        #
        # Launch the task that is named by
        #
        if 'factors' in cargs:
            config_list = factorial_iterator(cargs['factors'], cargs, config_parameters)
        else:
            config_list = factorial_iterator({}, cargs, config_parameters)

        for CONFIG in config_list:
            workflow = CONFIG['workflow']
            if workflow not in tasks:
                print("WARNING: Unknown workflow '%s' specified in block" % workflow)
                continue
            run_workflow(tasks[workflow], CONFIG, workflow_blocks, global_tempdir)


def run_block_parallel_factors(workflow_blocks, args, config_parameters):
    #
    # Serially execute all workflows, but parallelize the execution of workflow factors
    #
    tasks = registered_tasks()
    for cargs in workflow_blocks[args.block]:
        workflow = cargs.get('workflow', None)
        if workflow is None:
            print("WARNING: No workflow specified in '%s' block" % args.block)
            continue
        if not 'verbose' in cargs:
            cargs['verbose'] = args.verbose
        #
        # Launch the task that is named by
        #
        if 'factors' in cargs:
            config_list = list(factorial_iterator(cargs['factors'], cargs, config_parameters))
        else:
            config_list = list(factorial_iterator({}, cargs, config_parameters))

        if len(config_list) >= args.min_parallel_workflows:
            # Parallel
            good_configs = []
            for CONFIG in config_list:
                if workflow not in tasks:
                    print("WARNING: Unknown workflow '%s' specified in block" % workflow)
                    continue
                good_configs.append( CONFIG )
            with joblib.Parallel(n_jobs=args.np, verbose=args.verbose*10) as parallel:
                parallel( joblib.delayed(run_workflow)(tasks[CONFIG['workflow']], CONFIG, workflow_blocks, global_tempdir) for CONFIG in good_configs)
        else:
            # Serial
            for CONFIG in config_list:
                workflow = CONFIG['workflow']
                if workflow not in tasks:
                    print("WARNING: Unknown workflow '%s' specified in block" % workflow)
                    continue
                run_workflow(tasks[workflow], CONFIG, workflow_blocks, global_tempdir)

def run_block_parallel_workflows(workflow_blocks, args, config_parameters):
    #
    # Parallelize execution of workflows and workflow factors
    #
    tasks = registered_tasks()
    workflows = []
    for cargs in workflow_blocks[args.block]:
        workflow = cargs.get('workflow', None)
        if workflow is None:
            print("WARNING: No workflow specified in '%s' block" % args.block)
            continue
        if not 'verbose' in cargs:
            cargs['verbose'] = args.verbose
        #
        # Launch the task that is named by
        #
        if 'factors' in cargs:
            config_list = list(factorial_iterator(cargs['factors'], cargs, config_parameters))
        else:
            config_list = factorial_iterator({}, cargs, config_parameters)

        for CONFIG in config_list:
            if workflow not in tasks:
                print("WARNING: Unknown workflow '%s' specified in block" % workflow)
                continue
            workflows.append( CONFIG )

        if len(workflows) >= args.min_parallel_workflows:
            # Parallel
            with joblib.Parallel(n_jobs=args.np, verbose=args.verbose*10) as parallel:
                parallel( joblib.delayed(run_workflow)(tasks[CONFIG['workflow']], CONFIG, workflow_blocks, global_tempdir) for CONFIG in workflows)
        else:
            # Serial
            for CONFIG in workflows:
                workflow = CONFIG['workflow']
                run_workflow(tasks[workflow], CONFIG, workflow_blocks, global_tempdir)


def run_block_multiworkflow(CONFIG, warnings, config):
    for block in CONFIG['blocks']:
        if 'num_processes' in CONFIG:
            np = int(CONFIG['num_processes'])
            if CONFIG.get('parallel_workflows',False):
                args = Options(block=block,
                                np=np,
                                parallelize_workflow=True,
                                verbose=CONFIG['verbose'],
                                min_parallel_workflows=10)
                run_block_parallel_workflows(config, args, get_config_parameters())
            else:
                args = Options(block=block,
                                np=np,
                                parallelize_factors=True,
                                verbose=CONFIG['verbose'],
                                min_parallel_workflows=10)
                run_block_parallel_factors(config, args, get_config_parameters())
        else:
            args = Options(block=block, verbose=CONFIG['verbose'])
            run_block_serial(config, args, get_config_parameters())

    save_metadata(CONFIG, warnings)


def run(args):
    #
    #
    # Get help information
    #
    if args.help_workflows:
        tasks = registered_tasks()
        print("")
        print("Available Workflows")
        print("-------------------")
        for t in sorted(tasks.keys()):
            print(t)
            print("  "+tasks[t].description)
        return
    #
    # Load the YAML configuration file
    #
    workflow_blocks = load_configfile(args.config_file)
    if args.verbose:
        print("Configuration Arguments")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(get_config_parameters())
        pp.pprint(workflow_blocks)
        print("")
    #
    # Iterate over all workflow blocks for the specified subcommand
    #
    if not (args.block in workflow_blocks):
        print("WARNING: No workflow block '%s' specified" % args.block)
        return

    if not joblib_available and (args.parallelize_workflows or args.parallelize_factors):
        print("ERROR: Cannot parallelize workflows.  The 'joblib' packages is not installed.")
        return

    if args.parallelize_workflows:
        run_block_parallel_workflows(workflow_blocks, args, get_config_parameters())

    elif args.parallelize_factors:
        run_block_parallel_factors(workflow_blocks, args, get_config_parameters())

    else:
        run_block_serial(workflow_blocks, args, get_config_parameters())

