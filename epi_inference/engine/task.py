__all__ = ['Task']


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

