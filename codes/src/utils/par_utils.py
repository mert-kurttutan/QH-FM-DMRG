class SimpleNamespace(object):
    '''
    Useful when defining a set of parameters
    '''
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)