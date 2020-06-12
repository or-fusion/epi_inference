__all__ = ['save_metadata', 'compute_basename']

import yaml
import datetime


def compute_basename(filename):
    if filename.endswith('.csv') or filename.endswith('.yml') or filename.endswith('jsn'):
        return filename[:-4]
    if filename.endswith('.yaml') or filename.endswith('.json'):
        return filename[:-5]
    return None


def save_metadata(cargs, warnings, data=None):
        #
        # If the YAML data contains an 'output*' file, then 
        # create a YAML file with metadata
        #
        ofname = None
        ofname = cargs.get('output',ofname)
        ofname = cargs.get('output_csv',ofname)
        ofname = cargs.get('output_json',ofname)
        if ofname is not None:
            metadata = {}
            metadata['timestamp'] = str(datetime.datetime.now())
            metadata['configuration'] = {}
            for key in cargs:
                metadata['configuration'][key] = cargs[key]
            if not data is None and len(data) > 0:
                metadata['workflow parameters'] = data
            if len(warnings) > 0:
                metadata['warnings'] = warnings
            #
            basename = compute_basename(ofname)
            if basename is None:
                raise RuntimeError("Unknown output suffix: "+ofname)
            dfile = basename+"_meta.yml"
            #
            print("Writing file: "+dfile)
            with open(dfile, 'w') as OUTPUT:
                yaml.dump(metadata, OUTPUT)

