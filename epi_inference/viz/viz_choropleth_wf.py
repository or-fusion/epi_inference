__all__ = ['run', 'inference_choropleth']

import json
import os.path

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..viz.choropleth import create_us_choropleth


def run(CONFIG, warnings):
    if not 'input_json' in CONFIG:
        print("No 'input_json' value specified")
        warnings.append("No 'input_json' value specified")
    elif not os.path.exists(CONFIG['input_json']):
        print("File %s does not exist" % CONFIG['input_json'])
        warnings.append("File %s does not exist" % CONFIG['input_json'])
    else:
        with open(CONFIG['input_json'], 'r') as INPUT:
            results_json = json.load(INPUT)
        create_us_choropleth(results_json=results_json,
            value_key='beta',
            output_html=CONFIG['output_html'],
            description=CONFIG.get('description','Choropleth Plot of Estimated COVID Transmission Rates'),
            show_browser=CONFIG.get('show_browser',None))


class Viz_Choropleth(Task):

    def __init__(self):
        Task.__init__(self, "viz_choropleth",
            "Create visualization of estimated beta values.")

    def validate(self, CONFIG):
        valid_options = set(['description', 'input_json', 'output_html', 'show_browser', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(Viz_Choropleth())

