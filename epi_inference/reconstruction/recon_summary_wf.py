__all__ = ['recon_summary']

import sys
import os.path
import json
import glob
import numpy as np
from pyutilib.misc import timing

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata


def summary_narrow(OUTPUT, input_json_files, scenario_index, counties):
    counties = set(counties)
    first = True
    series = ['transmissions', 'S', 'E', 'I1', 'I2', 'I3', 'R']
    values = []
    data = {}

    for filename in glob.glob(input_json_files):
        if not os.path.exists(filename):
            raise RuntimeError("ERROR: Reconstruction JSON file does not exist: "+ filename)
        #
        with open(filename,'r') as INPUT:
            raw = json.load(INPUT)
        if first:
            #
            # Process the first JSON file
            #
            if len(counties) == 0:
                counties = list(sorted(raw.keys()))
            for fips in raw:
                curr = raw[fips]
                if 'E' not in curr:
                    continue
                for key in sorted(curr.keys()):
                    if key in series or key in values or key == 'dates':
                        continue
                    elif key != "FIPS":
                        values.append( key )
                break
            #
            values.remove(scenario_index)
            OUTPUT.write("fips,"+",".join(values)+",date,series,mean,Q25,Q50,Q75")
            OUTPUT.write("\n")
            first=False

        for fips in counties:
            if not fips in raw:
                continue
            _value = tuple(raw[fips][val] for val in values)
            for d in range(len(raw[fips]['dates'])):
                for s in series:
                    dateval = raw[fips]['dates'][d]
                    if (fips,_value,dateval) not in data:
                        data[fips,_value,dateval] = {}
                    if s not in data[fips,_value,dateval]:
                        data[fips,_value,dateval][s] = []
                    data[fips,_value,dateval][s].append( float(raw[fips][s][d]) )

        sys.stdout.write(".")
    sys.stdout.write("\n")

    for key in sorted(data.keys()):
        fips, _value, d = key

        for s in series:
            vals = data[key][s]
            mean = np.mean(vals)
            quartiles = np.quantile(vals, [0.25, 0.5, 0.75])

            prefix = list(map(str,[fips]+list(_value)+[d,s]))
            results = [str(mean), str(quartiles[0]), str(quartiles[1]), str(quartiles[2])]
            OUTPUT.write(",".join(prefix + results))
            OUTPUT.write("\n")


def recon_summary(input_json_files, output_csv, scenario_index, counties=None):
    # Write CSV file
    print("Writing reconstruction summary: "+output_csv)

    with open(output_csv,'w') as OUTPUT:
        summary_narrow(OUTPUT, input_json_files, scenario_index, counties)


class Recon_Summary_Workflow(Task):

    def __init__(self):
        Task.__init__(self, "recon_summary",
            "Summarize reconstructions in a CSV file.")

    def validate(self, CONFIG):
        valid_options = set(['scenario_index', 'input_json', 'output_csv', 'counties', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        recon_summary(
                CONFIG['input_json'],
                counties=CONFIG.get('counties',[]),
                output_csv=CONFIG['output_csv'],
                scenario_index=CONFIG['scenario_index'])


register_task(Recon_Summary_Workflow())

