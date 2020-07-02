__all__ = ['run', 'inference_json2csv']

import sys
import os.path
try:
    import ujson as json
except:
    import json
import csv

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata


def inference_json2csv(input_json, output_csv=None, datadir=None):
    if datadir:
        full_infile = os.path.join(datadir, input_json)
    else:
        full_infile = input_json
    if not os.path.exists(full_infile):
        raise RuntimeError("ERROR: Inference JSON file does not exist: "+ full_infile)
    #
    print("Processing inference file: "+full_infile)
    with open(full_infile,'r') as INPUT:
        raw = json.load(INPUT)

    # Figure out CSV output filename
    if output_csv:
        pass
    elif input_json.endswith('jsn'):
        output_csv = input_json[:-4]+".csv"
    elif input_json.endswith('json'):
        output_csv = input_json[:-5]+".csv"
    else:
        raise RuntimeError("ERROR: Cannot infer CSV output file name")
    if datadir:
        full_outfile = os.path.join(datadir,output_csv)
    else:
        full_outfile = output_csv

    # Write CSV file
    print("Writing results summary: "+full_outfile)
    with open(full_outfile,'w') as OUTPUT:
        sorted_fips = list(sorted(raw.keys()))
        OUTPUT.write("FIPS,Date,Beta,Status\n")
        for fips in sorted_fips:
            for d in range(len(raw[fips]['date'])):
                if raw[fips]['beta'][d] is None:
                    OUTPUT.write('"%s",%s,,%s\n' % (fips, raw[fips]['date'][d], raw[fips]['status'][d]))
                else:
                    OUTPUT.write('"%s",%s,%f,%s\n' % (fips, raw[fips]['date'][d], raw[fips]['beta'][d], raw[fips]['status'][d]))


def run(CONFIG, warnings):
    inference_json2csv(CONFIG['input_json'], output_csv=CONFIG.get('output_csv',None), datadir=CONFIG.get('datadir',None))


class Inference_JSON2CSV_Workflow(Task):

    def __init__(self):
        Task.__init__(self, "inference_json2csv",
            "Convert inference JSON to CSV.")

    def run(self, data, CONFIG):
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(Inference_JSON2CSV_Workflow())

