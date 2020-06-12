__all__ = ['run']

import sys
import pandas as pd

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata

from .load_results import load_df_casedata
from .misc import save_collect


def run(CONFIG, warnings):
    #
    # Process county case data.  We assume that data is organized within
    # a single directory, where each CSV file reports case data for a single
    # county.  See the epi_inference/examples/countydata directory structure for
    # an example.
    #
    df = load_df_casedata(CONFIG["files"],
                                datadir=CONFIG["dir"],
                                datacol=CONFIG.get("column", None),
                                days_before_first=CONFIG.get("days_before_first", None),
                                days_after_first=CONFIG.get("days_after_first",None))

    if df is None:                                  # pragma: no cover
        print("ERROR: no case data loaded")
        sys.exit(1)
    if CONFIG['verbose']:
        print("Data Summary")
        print(df)
        print("")

    save_collect(CONFIG['output'], df, CONFIG['verbose'], warnings)
    save_metadata(CONFIG, warnings)


class CollectCaseData(Task):

    def __init__(self):
        Task.__init__(self, "casedata",
            "Collect case data into a single CSV file for reconstruction and inference.")

    def validate(self, CONFIG):
        valid_options = set(['files', 'dir', 'column', 'days_before_first', 'days_after_first', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        run(CONFIG, self._warnings)


register_task(CollectCaseData())

