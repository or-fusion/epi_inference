import os
#import itertools
#from pyutilib.misc import Options
#import string
import datetime
import json
import pandas as pd
from deepdiff import DeepDiff


def roundall(*args):
    """
    Pass in any number of lists as arguments, and this function
    will return new lists (in the same order) with all values
    rounded to integers
    """
    ret_lists = list()
    for arg in args:
        assert type(arg) is list
        l = [None]*len(arg)
        for i,v in enumerate(arg):
            l[i] = round(v)
        ret_lists.append(l)
    return ret_lists


#
# NOTE: This assumes that it's safe to ignore the time when encoding a datetime object
#
class ToStr_JSONEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return obj.__str__()
        elif isinstance(obj, datetime.date):
            return obj.__str__()
        #
        # Convert numpy int/float types
        #
        try:
            return int(obj)
        except:
            try:
                return float(obj)
            except:
                pass
        return json.JSONEncoder.default(self, obj)

    def _encode(self, obj):
        def transform_date(o):
            return self._encode(o.strftime("%Y-%m-%d") if isinstance(o, datetime.datetime) or isinstance(o, datetime.date) else o)
        if isinstance(obj, dict):
            return {transform_date(k): transform_date(v) for k, v in obj.items()}
        elif isinstance(obj, list) or isinstance(obj, set):
            return [transform_date(l) for l in obj]
        else:
            return obj

    def encode(self, obj):
        return super(ToStr_JSONEncoder, self).encode(self._encode(obj))


def load_population(input_csv, index):
    try:
        population_df = pd.read_csv(input_csv, encoding="ISO-8859-1", dtype={index:'str'})
        population_df = population_df.set_index(index)
    except:                                                         # pragma: no cover
        raise RuntimeError("ERROR reading file "+input_csv)
    return population_df


def save_results(results, output):
    print("Writing results in file "+output)
    filedir = os.path.dirname(output)
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    with open(output,'w') as OUTPUT:
        json.dump(results, OUTPUT, cls=ToStr_JSONEncoder, indent=4)


def compare_csv(output, gold, index_col=None, check_exact=False, sort=True):
    if index_col is None:
        outputdf = pd.read_csv(output)
        golddf = pd.read_csv(gold)
    else:
        outputdf = pd.read_csv(output, index_col=index_col)
        golddf = pd.read_csv(gold, index_col=index_col)

    # the dataframes may be the same, but just in a different order
    if sort:
        columns = list(outputdf.columns)
        outputdf.sort_values(by=columns, inplace=True, ignore_index=True)
        golddf.sort_values(by=columns, inplace=True, ignore_index=True)
    pd.testing.assert_frame_equal(left=outputdf, right=golddf, check_exact=check_exact)
    return outputdf, golddf


def compare_json(output_file, baseline_file, significant_digits=8):            # pragma: no cover
    with open(output_file,'r') as INPUT:
        output = json.load(INPUT)
    with open(baseline_file,'r') as INPUT:
        baseline = json.load(INPUT)
    d = DeepDiff(baseline, output, significant_digits=significant_digits)
    if len(d) != 0:
        print('DIFFERENCES IN JSON')
        print(d)
    assert(len(d) == 0)
    return output, baseline

