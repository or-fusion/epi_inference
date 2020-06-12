__all__ = ['factorial_iterator', 'load_configfile']

import sys
import itertools
import string
import yaml
from pyutilib.misc import Options

from .config_parameters import get_config_parameters, set_config_parameters


def _safe_eval(value, levels):
    tmp = string.Template(value).safe_substitute(**levels)
    try:
        return int(tmp)
    except:
        try:
            return float(tmp)
        except:
            pass
    return tmp

def _eval(value, levels, other_values):
    tmp = string.Template(value).safe_substitute(**levels)
    tmp = string.Template(tmp).substitute(**other_values)
    try:
        return int(tmp)
    except:
        try:
            return float(tmp)
        except:
            pass
    return tmp


def factorial_iterator(factors, config, other_values={}):
    assert(config.get('factor_levels',None) is None)

    factor_names = list(sorted(factors.keys()))
    factor_list = list(factors[f] for f in factor_names)
    for fac_ in itertools.product(*factor_list):
        levels = {}
        CONFIG = Options()
        for i in range(len(fac_)):
            levels[factor_names[i]] = str(fac_[i])
            #CONFIG[factor_names[i]] = fac_[i]
        for key, value in config.items():
            if type(value) is str and '${' in value:
                CONFIG[key] = _eval(value, levels, other_values)
            elif type(value) is dict:
                CONFIG[key] = Options()
                for _key, _value in value.items():
                    if type(_value) is str and '${' in _value:
                        CONFIG[key][_key] = _eval(_value, levels, other_values)
                    else:
                        CONFIG[key][_key] = _value
            else:
                CONFIG[key] = value

        CONFIG.factors = None
        CONFIG.factor_levels = levels
        #
        # Yield a tuple with a modified configuration object
        #
        yield CONFIG


def yaml_parser__load_yaml_data(loader, node):
    value = loader.construct_scalar(node)
    value = string.Template(value).substitute(**get_config_parameters())
    with open(value, 'r') as INPUT:
        return yaml.load(INPUT, Loader=yaml.Loader)


def load_configfile(config_file):
    yaml.add_constructor(u'!LoadYAMLFile', yaml_parser__load_yaml_data)
    set_config_parameters({})

    with open(config_file, 'r') as INPUT:
        try:
            first = True
            configs = []
            for config in yaml.load_all(INPUT, Loader=yaml.Loader):
                if first:
                    first = False
                    for key in config:
                        if type(config[key]) in [list, dict]:
                            configs.append(config)
                            break
                    if len(configs) == 1:
                        continue
                    #
                    # It looks like this is a config block
                    # process is *now*, before the workflow block
                    # is processed.
                    #
                    # NOTE: This iterates at most 10 times trying to
                    # replace strings, and then gives up.
                    #
                    parameters = config
                    counter=10
                    while counter > 0:
                        flag=True
                        tmp = {}
                        for key,value in parameters.items():
                            if type(value) is str and '${' in value:
                                tmp[key] = _safe_eval(value, parameters)
                                flag = False
                            else:
                                tmp[key] = value
                        parameters = tmp
                        if flag:
                            break
                        counter = counter - 1
                    set_config_parameters(parameters)
                else:
                    configs.append(config)
        except yaml.YAMLError as exc:                   # pragma: nocover
            print("ERROR: problem parsing YAML file '%s'" % config_file)
            print(exc)
            sys.exit(1)
    #
    # Error checks
    #
    if len(configs) != 1:
        raise RuntimeError("Problem loading configuration file.  Can only process a YAML file with one or two YAML documents.")
    return configs[0]

