__all__ = ['register_task', 'registered_tasks']

#
# mapping from TaskName -> Task object instance
#
global_tasks = {}

def register_task(task):
    if task.name in global_tasks:
        raise RuntimeError("Task '%s' has already been registered!" % task.name)
    global_tasks[task.name] = task


def registered_tasks():
    return global_tasks
