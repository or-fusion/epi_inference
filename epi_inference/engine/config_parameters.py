__all__ = ['set_config_parameters', 'get_config_parameters']

global_config_parameters = {}

def set_config_parameters(config_parameters):
    global global_config_parameters
    global_config_parameters = config_parameters

def get_config_parameters():
    return global_config_parameters
