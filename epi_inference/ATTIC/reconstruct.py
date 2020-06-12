import csv
import json
import sys
#import argparse
import yaml
import os.path
#import shutil
import datetime
import pprint
from .util import factorial_iterator, ToStr_JSONEncoder
from .formulations import reconstruction as recon
from pyutilib.misc import Options
import pandas as pd


def process_config(cfg):
    config = Options()

    #config.formulation = cfg['formulation']
    #config.ntrials = cfg.get('ntrials', None)
    config.input_csv = cfg['input_csv']
    config.population = cfg.get('population', None)
    config.county = cfg.get('county', None)
    #config.column = cfg.get('column', None)
    config.reporting_factor = cfg.get('reporting_factor', 1.0)
    config.deltaP = cfg.get('deltaP', 7)
    config.sigma = cfg.get('sigma', None)
    config.gamma = cfg.get('gamma', None)
    #config.factor_levels = cfg.get('factor_levels', None)
    #config.bootstrap = cfg.get('bootstrap', Options())
    #config.analysis_window = cfg.get('analysis_window', Options())

    # TODO - deprecate the use of the geodata CSV file option
    # TODO - deprecate the use of the geodata CSV file option
    if 'population_csv' not in cfg:
        config.population_csvfile = cfg.get('geodata_csv', config.input_csv[:-4] + "_geodata.csv")
        config.population_csvcolumn = 'pop2010'
        config.population_csvindex = 'geoid'
    else:
        config.population_csvfile = cfg['population_csv']['file']
        config.population_csvcolumn =  cfg['population_csv']['population']
        config.population_csvindex =  cfg['population_csv']['index']

    return config

def run_reconstruction(df, population_df, CONFIG, verbose):
    population_config = CONFIG.population
    population_csvcolumn = CONFIG.population_csvcolumn
    sigma = CONFIG.sigma
    gamma = CONFIG.gamma
    deltaP = CONFIG.deltaP
    county = CONFIG.county
    reporting_factor = CONFIG.reporting_factor

    all_results = list()

    for t in df:
        if county is not None and t != county:
            continue

        if population_config is None:
            if t not in population_df[population_csvcolumn]:                # pragma: no cover
                print("WARNING: county "+str(t)+" does not have population data available.")
                continue
            population = population_df[population_csvcolumn][t]
        else:
            population = population_config
        cumulative_reported_cases = df[t].to_list()

        if df[t][-1] == 0:
            results = {}

        else:
            #
            # DO SIMULATION HERE
            #
            Cdates = pd.date_range(end=datetime.datetime(year=2020, month=4, day=12),
                                    periods=len(cumulative_reported_cases)).to_pydatetime().tolist()

            rdates, rcases, dates, T, S, E, I1, I2, I3, R = recon.reconstruct_states_deterministic_decay(
                                Cdates=Cdates,
                                cumulative_reported_cases=cumulative_reported_cases,
                                population=population,
                                sigma=sigma,
                                gamma=gamma/3,
                                reporting_factor=reporting_factor,
                                report_delay=deltaP
                                )

            results = {}
            results['rdates'] = rdates
            results['rcases'] = rcases
            results['dates'] = dates
            results['T'] = T
            results['S'] = S
            results['E'] = E
            results['I1'] = I1
            results['I2'] = I2
            results['I3'] = I3
            results['R'] = R
            results['population'] = population

        results['FIPS'] = t
        #
        # Collect results in a list
        #
        all_results.append( results )
    return all_results

def check_config(config):
    if config.county is not None and type(config.county) is not str:        # pragma: no cover
        print("ERROR: county id must be specified as a string")
        sys.exit(1)
    try:
        assert(os.path.exists(config.input_csv))
    except:                                                                 # pragma: no cover
        print("ERROR: input file "+config.input_csv+" does not exist")
        raise
    assert type(config.reporting_factor) is float


def run(args):
    try:
        with open(args.config_file, 'r') as INPUT:
            config = yaml.safe_load(INPUT)
    except yaml.YAMLError as exc:                   # pragma: nocover
        print("ERROR: problem parsing YAML file")
        print(exc)
        sys.exit(1)

    if 'reconstruct' not in config:
        raise ValueError('No "reconstruct" key found in the YAML config')

    for cfg in config.get('reconstruct', []):
        #
        # Process county case data, and execute reconstructions to predict the 
        # temporal evolution of the compartments in the epi model.
        #
        verbose = cfg.get('verbose', args.verbose)
        factors = cfg.get('factors', None)
        output = cfg.get('output_json', None)

        assert output is not None
        assert type(verbose) is bool

        config = process_config(cfg)
        if verbose:
            print("Configuration Arguments")
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(config)
            print("")

        all_results = []

        if factors is None:
            config_list = [config]
        else:
            config_list = factorial_iterator(factors, config)

        for CONFIG in config_list:
            try:
                population_df = pd.read_csv(CONFIG.population_csvfile, encoding="ISO-8859-1", dtype={CONFIG.                population_csvindex:'str'})
                population_df = population_df.set_index(CONFIG.population_csvindex)
            except:                                                         # pragma: no cover
                print("ERROR reading file "+CONFIG.population_csvfile)
                raise
            check_config(CONFIG)
            if CONFIG.population is None and CONFIG.county is not None:
                CONFIG.population = population_df[CONFIG.population_csvcolumn][CONFIG.county]
                if verbose:
                    print("County: %s    Population: %s" % (str(CONFIG.county), str(CONFIG.population)))

            #
            # Load the dataframe and experimental metadata (if it's available)
            #
            print("Input File: "+CONFIG.input_csv)
            df = pd.read_csv(CONFIG.input_csv, index_col='Date')
            #
            # Execute inference
            #
            results = run_reconstruction(df, population_df, CONFIG, verbose)
            #
            # Augment reported results
            #
            for trial in results:
                if CONFIG.factor_levels is not None:
                    for key, value in CONFIG.factor_levels.items():
                        if not key in trial:
                            trial[key] = value
                all_results.append( trial )

        #
        # Save results
        #
        print("Writing results in file "+output)
        filedir = os.path.dirname(output)
        if not os.path.exists(filedir):
            os.makedirs(filedir)
        with open(output,'w') as OUTPUT:
            json.dump(all_results, OUTPUT, cls=ToStr_JSONEncoder, indent=4)
        #
        # Create a YAML file with metadata
        #
        metadata = {}
        metadata['timestamp'] = str(datetime.datetime.now())
        metadata['configuration'] = {}
        for key in cfg:
            metadata['configuration'][key] = cfg[key]
        if output.endswith(".jsn"):
            metaoutput = output[:-4]+"_meta.yml"
        elif output.endswith(".json"):
            metaoutput = output[:-5]+"_meta.yml"
        else:
            metaoutput = output+"_meta.yml"
        print("Writing results metadata in file "+metaoutput)
        with open(metaoutput, 'w') as OUTPUT:
            yaml.dump(metadata, OUTPUT)
        #
        # Print data
        #
        #for r in all_results:
        #    print(r)

