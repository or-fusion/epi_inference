__all__ = ['run', 'inference_choropleth']

import json
import os.path

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..viz.choropleth import create_us_choropleth


def run(CONFIG, warnings):
    if not 'input_json' in CONFIG:
        warnings.append("No 'input_json' value specified")
    elif not os.path.exists(CONFIG['input_json']):
        warnings.append("FIle %s does not exist" % CONFIG['input_json'])
    else:
        with open(CONFIG['input_json'], 'r') as INPUT:
            results_json = json.load(INPUT)
        create_us_choropleth(results_json=results_json, value_key='beta')


class Viz_Choropleth(Task):

    def __init__(self):
        Task.__init__(self, "viz_choropleth",
            "Create visualization of estimated beta values.")

    def run(self, data, CONFIG):
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(Viz_Choropleth())

