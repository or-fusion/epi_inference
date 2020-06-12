__all__ = ['run']

import sys
import json
from pyutilib.misc import timing

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from ..util import load_population, save_results
from ..formulations.multinode_mobility_window_decay_lsq import run_multinode_mobility_window_decay_lsq
from ..formulations.multinode_mobility_window_decay_lsq_old import run_multinode_mobility_window_decay_lsq_old
from ..formulations.multinode_mobility_window_decay_lsq_poek import run_multinode_mobility_window_decay_lsq_poek


def run(CONFIG, warnings):
    #
    # Load the reconstruction data 
    #
    with open(CONFIG['reconstruction_json'],'r') as INPUT:
        recon = json.load(INPUT)
    #
    # Load the mobility data 
    #
    with open(CONFIG['mobility_json'],'r') as INPUT:
        mobility = json.load(INPUT)
    #
    # Perform inference
    #
    #try:
    if True:
        if CONFIG.get('version','new') == 'new':
            results = run_multinode_mobility_window_decay_lsq(recon=recon, mobility=mobility, analysis_window=CONFIG['analysis_window'], verbose=CONFIG['verbose'])
        elif CONFIG.get('version','new') == 'poek':
            results = run_multinode_mobility_window_decay_lsq_poek(recon=recon, mobility=mobility, analysis_window=CONFIG['analysis_window'], verbose=CONFIG['verbose'])
        else:
            results = run_multinode_mobility_window_decay_lsq_old(recon=recon, mobility=mobility, analysis_window=CONFIG['analysis_window'], verbose=CONFIG['verbose'])
    #except Exception as err:
    else:
        print("ERROR: Unexpected exception '%s'" % str(err))
        results = {}
        warnings.append(str(err))
    #
    # Save results
    #
    save_results(results, CONFIG['output_json'])
    save_metadata(CONFIG, warnings)


class InferenceMobilityWindows(Task):

    def __init__(self):
        Task.__init__(self, "estimate_beta_windows_with_mobility",
            "Estimate beta over different time windows using inter-county mobility information.")

    def validate(self, args):
        pass

    def run(self, data, CONFIG):
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(InferenceMobilityWindows())

