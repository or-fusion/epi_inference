__all__ = ['run']

import sys
try:
    import rpy2
    rpy2_available = True
except:
    rpy2_available = False

from ..engine.task import Task
from ..engine.task_registry import register_task
from ..engine.misc import save_metadata


def run(CONFIG, warnings):
    print("Running test.R")
    import rpy2.robjects as robjects
    robjects.r.source('../../../R_utilities/test.R')


class RTestWorkflow(Task):

    def __init__(self):
        Task.__init__(self, "Rtest",
            "Simple test of R interface from Python.")

    def run(self, data, CONFIG):
        if not rpy2_available:
            raise RuntimeError("ERROR: cannot execute Rtest workflow.  Package 'rpy2' is not available.")
        self._warnings = []
        run(CONFIG, self._warnings)

    def warnings(self):
        return self._warnings


register_task(RTestWorkflow())

