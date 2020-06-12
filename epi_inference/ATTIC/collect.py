import csv
import sys
import argparse
import yaml
import os.path
import shutil
import datetime
import pprint
from .util import factorial_iterator

from .load_results import load_df_expdata, load_df_casedata


def save_output(cargs, df, verbose, data=None):
    if verbose:
        print("Data Summary")
        print(df)
        print("")
   
    ofname = cargs["output"]
    print("Writing file: "+ofname)
    filedir = os.path.dirname(ofname)
    if filedir and not os.path.exists(filedir):         # pragma: no cover
        os.makedirs(filedir)
    if verbose and os.path.exists(ofname):              # pragma: no cover
        print("  WARNING: over-writing file "+ofname)
    df.to_csv(ofname, quoting=csv.QUOTE_NONNUMERIC)
    #
    if os.path.exists(os.path.join(cargs["dir"], "geodata.csv")):
        shutil.copyfile(os.path.join(cargs["dir"], "geodata.csv"), ofname[:-4]+"_geodata.csv")        
    #
    # Create a YAML file with metadata
    #
    metadata = {}
    metadata['timestamp'] = str(datetime.datetime.now())
    metadata['configuration'] = {}
    for key in cargs:
        metadata['configuration'][key] = cargs[key]
    if not data is None and len(data) > 0:
        metadata['simulation parameters'] = data
    dfile = ofname[:-4]+"_meta.yml"
    print("Writing file: "+dfile)
    with open(dfile, 'w') as OUTPUT:
        yaml.dump(metadata, OUTPUT)

def run(args):
    #
    # Load the YAML configuration file
    #
    with open(args.config_file, 'r') as INPUT:
        try:
            config = yaml.safe_load(INPUT)
        except yaml.YAMLError as exc:                   # pragma: nocover
            print("ERROR: problem parsing YAML file")
            print(exc)
            sys.exit(1)
    if args.verbose:
        print("Configuration Arguments")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(config)
        print("")

    if 'expdata' in config:
        #
        # Process experimental data.  We assume that data is organized
        # within a subdirectory in separate directories named 'exp<id>'.
        # See the epi_inference/examples/expdata directory structure for
        # an example.
        #
        cargs = config['expdata']
        for i in range(len(cargs)):
            if 'factors' in cargs[i]:
                config_list = factorial_iterator(cargs[i]['factors'], cargs[i])
            else:
                config_list = [cargs[i]]

            for CONFIG in config_list:
                df, data = load_df_expdata(expdir=CONFIG["dir"],
                                                county=CONFIG["county"],
                                                trial=CONFIG.get('trial',None),
                                                days_before_first=CONFIG.get("days_before_first", None),
                                                days_after_first=CONFIG.get("days_after_first",None),
                                                expnum=CONFIG["expnum"])

                if df is None:                                      # pragma: no cover
                    print("ERROR: no experimental data loaded")
                    sys.exit(1)
                save_output(CONFIG, df, args.verbose, data=data)

    elif 'casedata' in config:
        #
        # Process county case data.  We assume that data is organized within
        # a single directory, where each CSV file reports case data for a single
        # county.  See the epi_inference/examples/countydata directory structure for
        # an example.
        #
        cargs = config['casedata']
        for i in range(len(cargs)):
            if 'factors' in cargs[i]:
                config_list = factorial_iterator(cargs[i]['factors'], cargs[i])
            else:
                config_list = [cargs[i]]

            for CONFIG in config_list:
                df = load_df_casedata(CONFIG["files"],
                                            datadir=CONFIG["dir"],
                                            datacol=CONFIG.get("column", None),
                                            days_before_first=CONFIG.get("days_before_first", None),
                                            days_after_first=CONFIG.get("days_after_first",None))

                if df is None:                                  # pragma: no cover
                    print("ERROR: no case data loaded")
                    sys.exit(1)
                save_output(CONFIG, df, args.verbose)

