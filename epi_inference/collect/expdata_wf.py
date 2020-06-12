__all__ = ['run']

import os
import pandas as pd
import shutil

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata, compute_basename

from .load_results import load_df_expdata
from .misc import save_collect


def run(CONFIG, warnings):
    #
    # Process county case data.  We assume that data is organized within
    # a single directory, where each CSV file reports case data for a single
    # county.  See the epi_inference/examples/countydata directory structure for
    # an example.
    #
    df, data = load_df_expdata(expdir=CONFIG["dir"],
                                    county=CONFIG["county"],
                                    trial=CONFIG.get('trial',None),
                                    days_before_first=CONFIG.get("days_before_first", None),
                                    days_after_first=CONFIG.get("days_after_first",None),
                                    expnum=CONFIG["expnum"])

    if df is None:                                      # pragma: no cover
        raise RuntimeError("ERROR: no experimental data loaded")
    if CONFIG['verbose']:
        print("Data Summary")
        print(df)
        print("")

    save_collect(CONFIG['output'], df, CONFIG['verbose'], warnings)
    save_metadata(CONFIG, warnings, data=data)

    if os.path.exists(os.path.join(CONFIG["dir"], "geodata.csv")):
        shutil.copyfile(os.path.join(CONFIG["dir"], "geodata.csv"), compute_basename(CONFIG['output'])+"_geodata.csv")



class CollectExpData(Task):

    def __init__(self):
        Task.__init__(self, "expdata",
            "Collect experiment data into a single CSV file for reconstruction and inference.")

    def validate(self, CONFIG):
        valid_options = set(['dir', 'county', 'trial', 'days_before_first', 'days_after_first', 'expnum', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run(CONFIG, self._warnings)


register_task(CollectExpData())

