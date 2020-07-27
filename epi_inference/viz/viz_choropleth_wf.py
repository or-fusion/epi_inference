__all__ = ['run_scenario', 'run_summary']

import os.path
import glob
import csv
try:
    import ujson as json
except:
    import json

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..viz.choropleth import create_us_choropleth_scenario, create_us_choropleth_summary


def run_scenario(CONFIG, warnings):
    if not 'input_json' in CONFIG:
        print("No 'input_json' value specified")
        warnings.append("No 'input_json' value specified")
    elif not os.path.exists(CONFIG['input_json']):
        print("File %s does not exist" % CONFIG['input_json'])
        warnings.append("File %s does not exist" % CONFIG['input_json'])
    else:
        with open(CONFIG['input_json'], 'r') as INPUT:
            results_json = json.load(INPUT)
        create_us_choropleth_scenario(results_json=results_json,
            value_key='beta',
            output_html=CONFIG['output_html'],
            description=CONFIG.get('description','Choropleth Plot of Estimated COVID Transmission Rates'),
            show_browser=CONFIG.get('show_browser',None))

def run_summary(CONFIG, warnings):
    if not 'input_csv' in CONFIG:
        print("No 'input_csv' value specified")
        warnings.append("No 'input_csv' value specified")
    else:
        data = {}
        for filename in glob.glob(CONFIG['input_csv']):
            with open(filename, 'r') as INPUT:
                first = None
                last = None
                for row in csv.reader(INPUT):
                    if first is None:
                        first = row
                    last = row
                FIPS = filename.split('.')[-2].split('_')[-1]
                # NOTE - If we allow the user to specify the date, then we need to keep all the rows...
                data[FIPS] = [{first[i]:last[i] for i in range(len(first))}]
        create_us_choropleth_summary(summary_csv=data,
            value_key='qmean_filtered_est_beta',
            output_html=CONFIG['output_html'],
            description=CONFIG.get('description','Choropleth Plot of Mean Values of Estimated COVID Transmission Rates'),
            show_browser=CONFIG.get('show_browser',None))


class Viz_ChoroplethScenario(Task):

    def __init__(self):
        Task.__init__(self, "viz_choropleth_scenario",
            "Create visualization of a single set of estimated beta values.")

    def validate(self, CONFIG):
        valid_options = set(['description', 'input_json', 'output_html', 'show_browser', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run_scenario(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


class Viz_ChoroplethSummary(Task):

    def __init__(self):
        Task.__init__(self, "viz_choropleth_summary",
            "Create visualization of a summary of estimated beta values over a set of scenarios.")

    def validate(self, CONFIG):
        valid_options = set(['description', 'input_csv', 'output_html', 'show_browser', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run_summary(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(Viz_ChoroplethScenario())
register_task(Viz_ChoroplethSummary())

