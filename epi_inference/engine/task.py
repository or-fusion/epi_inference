__all__ = ['Task', 'GenericWorkflow']


class Task(object):

    def __init__(self, name, description):
        self.name = name
        self.description = description

    def validate(self, args):
        pass

    def run(self, data, args):
        pass

    def warnings(self):
        pass


class GenericWorkflow(Task):

    def __init__(self, name, description):
        Task.__init__(name, description)
        self.tasks = []

    def add_task(self, task):
        self.tasks.append(task)

    def validate(self, args):
        self.validate_workflow(args)
        for task in self.tasks:
            validate.validate(args)

    def run(self, data, args):
        for task in self.tasks:
            data = task.run(data, args)

    def validate_workflow(self, args):
        pass
