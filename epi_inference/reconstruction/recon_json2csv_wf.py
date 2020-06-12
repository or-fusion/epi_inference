__all__ = ['recon_single_json2csv', 'recon_many_json2csv']

import sys
import os.path
import json
import csv
import glob

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata


def write_wide(OUTPUT, raw, counties):
    counties = set(counties)
    OUTPUT.write('comp,')
    sorted_fips = list(sorted(raw.keys()))
    for fips in sorted_fips:
        if fips not in counties:
            continue
        OUTPUT.write('"%s",' % fips)
    OUTPUT.write("time\n")
    for d in range(len(raw[fips]['dates'])):
        for s in ['S', 'E', 'I1', 'I2', 'I3', 'R']:
            OUTPUT.write('"%s",' % s)
            for fips in sorted_fips:
                if fips not in counties:
                    continue
                OUTPUT.write("%s," % str(raw[fips][s][d]))
            OUTPUT.write('"%s"\n' % raw[fips]['dates'][d])


def write_flattened(OUTPUT, input_json_files, counties):
    counties = set(counties)
    first = True
    series = ['dates', 'transmissions', 'S', 'E', 'I1', 'I2', 'I3', 'R']
    values = []

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
                    if key in series or key in values:
                        continue
                    elif key != "FIPS":
                        values.append( key )
                break
            #
            OUTPUT.write("FIPS,"+",".join(values+series))
            OUTPUT.write("\n")
            first=False

        rows = []
        for fips in counties:
            if not fips in raw:
                continue
            for d in range(len(raw[fips]['dates'])):
                row = []
                row.append('"%s"' % fips)
                for val in values:
                    row.append('%s' % str(raw[fips][val]))
                for s in series:
                    row.append('%s' % str(raw[fips][s][d]))
                rows.append(",".join(row))
        OUTPUT.write("\n".join(rows))
        OUTPUT.write("\n")
        sys.stdout.write(".")
        sys.stdout.flush()

    sys.stdout.write("\n")


def write_narrow(OUTPUT, input_json_files, counties):
    counties = set(counties)
    first = True
    series = ['transmissions', 'S', 'E', 'I1', 'I2', 'I3', 'R']
    values = []

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
            OUTPUT.write("fips,"+",".join(values)+",date,series,value")
            OUTPUT.write("\n")
            first=False

        rows = []
        for fips in counties:
            if not fips in raw:
                continue
            prefix = []
            prefix.append('"%s"' % fips)
            for val in values:
                prefix.append('%s' % str(raw[fips][val]))
            for d in range(len(raw[fips]['dates'])):
                for s in series:
                    appendix = ['%s' % str(raw[fips]['dates'][d]), 
                                s,
                                '%s' % str(raw[fips][s][d])]
                    rows.append(",".join(prefix+appendix))
        OUTPUT.write("\n".join(rows))
        OUTPUT.write("\n")
        sys.stdout.write(".")
        sys.stdout.flush()

    sys.stdout.write("\n")


def recon_single_json2csv(input_json, output_csv=None, datadir=None, csv_format=None, counties=None):
    if datadir:
        full_infile = os.path.join(datadir, input_json)
    else:
        full_infile = input_json
    if not os.path.exists(full_infile):
        raise RuntimeError("ERROR: Reconstruction JSON file does not exist: "+ full_infile)
    #
    print("Processing reconstruction file: "+full_infile)
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
    if csv_format == 'wide':
        with open(full_outfile,'w') as OUTPUT:
            write_wide(OUTPUT, raw, counties)
    elif csv_format == 'flatten':
        with open(full_outfile,'w') as OUTPUT:
            write_flattened(OUTPUT, full_infile, counties)
    else:
        with open(full_outfile,'w') as OUTPUT:
            write_narrow(OUTPUT, full_infile, counties)


def recon_many_json2csv(input_json_files, output_csv=None, datadir=None, counties=None, csv_format=None):
    # Figure out CSV output filename
    if datadir:
        full_outfile = os.path.join(datadir,output_csv)
    else:
        full_outfile = output_csv

    # Write CSV file
    print("Writing results summary: "+full_outfile)
    if csv_format == 'flatten':
        with open(full_outfile,'w') as OUTPUT:
            write_flattened(OUTPUT, input_json_files, counties)
    else:
        with open(full_outfile,'w') as OUTPUT:
            write_narrow(OUTPUT, input_json_files, counties)


class Recon_JSON2CSV_Workflow(Task):

    def __init__(self):
        Task.__init__(self, "recon_json2csv",
            "Convert reconstruction JSON to CSV.")

    def validate(self, CONFIG):
        valid_options = set(['format', 'input_json', 'counties', 'output_csv', 'datadir', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        recon_single_json2csv(
                CONFIG['input_json'],
                output_csv=CONFIG.get('output_csv',None),
                datadir=CONFIG.get('datadir',None),
                counties=CONFIG.get('counties',[]),
                csv_format=CONFIG.get('format','narrow'))


register_task(Recon_JSON2CSV_Workflow())


class Recon_Many_JSON2CSV_Workflow(Task):

    def __init__(self):
        Task.__init__(self, "recon_many_json2csv",
            "Create a CSV file summarizing reconstructions in JSON files.")

    def validate(self, CONFIG):
        valid_options = set(['format', 'input_json', 'counties', 'output_csv', 'datadir', 'verbose', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        recon_many_json2csv(
                CONFIG['input_json'],
                counties=CONFIG.get('counties',[]),
                output_csv=CONFIG.get('output_csv',None),
                datadir=CONFIG.get('datadir',None),
                csv_format=CONFIG.get('format','narrow'))


register_task(Recon_Many_JSON2CSV_Workflow())

