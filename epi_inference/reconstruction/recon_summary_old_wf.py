__all__ = ['recon_summary']

import sys
import os.path
try:
    import ujson as json
except:
    import json
import csv
import glob
import numpy as np
import pandas as pd
from pyutilib.misc import timing

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata


def summary_narrow(OUTPUT, reader, scenario_index):
    header = next(reader)
    offset = {header[i]:i for i in range(len(header))}
    #print("HEADER:",header)
    if not scenario_index in header:
        raise RuntimeError("The scenario index is not specified in the CSV file: '%s'" % scenario_index)
    header.remove('value')
    header.remove(scenario_index)
    #
    # Collect data from dataframe
    #
    timing.tic()
    data = {}
    for row in reader:
        key = tuple(row[offset[name]] for name in header)
        if key not in data:
            data[key] = []
        data[key].append( float(row[offset['value']]) )
    timing.toc(".  Collected data")
    #
    # Create summary CSV
    #
    OUTPUT.write(",".join(header+['mean','Q25', 'Q50', 'Q75']))
    OUTPUT.write("\n")
    for key, vals in data.items():
        #print(key,vals)
        #mean = statistics.mean(vals)
        mean = np.mean(vals)
        #quartiles = statistics.quantiles(vals, method='inclusive')
        quartiles = np.quantile(vals, [0.25, 0.5, 0.75])
        
        results = (str(mean), str(quartiles[0]), str(quartiles[1]), str(quartiles[2]))
        OUTPUT.write(",".join(map(str,key + results)))
        OUTPUT.write("\n")
    timing.toc(".  Wrote summary")


def recon_summary(input_csv, output_csv, scenario_index, csv_format=None):
    if not os.path.exists(input_csv):
        raise RuntimeError("ERROR: Reconstruction CSV file does not exist: "+ input_csv)
    #
    #print("Processing reconstruction file: "+input_csv)
    #with open(input_csv,'r') as INPUT:
        #df = pd.read_csv(INPUT)

    # Write CSV file
    print("Writing results summary: "+output_csv)
    if csv_format == 'wide':
        #with open(full_outfile,'w') as OUTPUT:
        #    write_wide(OUTPUT, raw, counties)
        pass
    elif csv_format == 'flatten':
        #with open(full_outfile,'w') as OUTPUT:
        #    write_flattened(OUTPUT, full_infile, counties)
        pass
    else:
        with open(output_csv,'w') as OUTPUT:
            with open(input_csv,'r') as INPUT:
                reader = csv.reader(INPUT)
                summary_narrow(OUTPUT, reader, scenario_index)


class Recon_Summary_WorkflowOLD(Task):

    def __init__(self):
        Task.__init__(self, "recon_summary_old",
            "Summarize a reconstruction CSV.")

    def validate(self, CONFIG):
        valid_options = set(['format', 'scenario_index', 'input_csv', 'output_csv', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        recon_summary(
                CONFIG['input_csv'],
                output_csv=CONFIG['output_csv'],
                scenario_index=CONFIG['scenario_index'],
                csv_format=CONFIG.get('format','narrow'))


register_task(Recon_Summary_WorkflowOLD())

