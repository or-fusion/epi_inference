__all__ = []

from . import driver
from .task import Task
from .task_registry import register_task


class BlockMultiWorkflow(Task):

    def __init__(self):
        Task.__init__(self, "execute_blocks",
            "Define a sequence of YAML blocks that are sequentially executed.  But parallelization options are applied to each in turn.")

    def validate(self, CONFIG):
        valid_options = set(['blocks', 'num_processes', 'parallel_workflows', 'verbose', 'output', 'factors', 'factor_levels', 'workflow'])
        for key in CONFIG:
            if key not in valid_options:
                raise RuntimeError("Unexpected configuration option: '%s'" % key)

    def run(self, data, CONFIG):
        self._warnings = []
        self.validate(CONFIG)
        driver.run_block_multiworkflow(CONFIG, self._warnings, data)


register_task(BlockMultiWorkflow())

