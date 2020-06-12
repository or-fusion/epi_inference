import argparse
import yaml
from .task_registry import registered_tasks


def driver():
    Parser = argparse.ArgumentParser(description='inference models')
    Parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Verbosity flag')
    subparsers = Parser.add_subparsers(title='subcommands', help="Help", dest="subparser_name") 
    #
    # Create subparsers
    #
    tasks = registered_tasks()
    parsers = []
    for name in sorted(tasks.keys()):
        parser = subparsers.add_parser(name, help=tasks[name].description)
        parser.add_argument('config_file', help='YAML configuration file')
        parsers.append(parser)
    #
    # Parse sys.argv
    #
    args = Parser.parse_args()
    if args.subparser_name is None:
        Parser.print_help()
        return
    #
    # Load the YAML configuration file
    #
    with open(args.config_file, 'r') as INPUT:
        try:
            config = yaml.safe_load(INPUT)
        except yaml.YAMLError as exc:                   # pragma: nocover
            print("ERROR: problem parsing YAML file")
            print(exc)
            sys.exit(1)
    if args.verbose:
        print("Configuration Arguments")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(config)
        print("")
    #
    # Iterate over all config blocks for the specified subcommand
    #
    if not (args.subparser_name in config):
        print("WARNING: No configuration blocks specified for '%s' subcommand" % args.subparser_name)

    for cargs in config[args.subparser_name]:
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
            #
            if ofname.endswith('.csv') or ofname.endswith('.jsn') or ofname.endswith('yml'):
                dfile = cargs['output'][:-4]+"_meta.yml"
            elif ofname.endswith('json') or ofname.endswith('yaml'):
                dfile = cargs['output'][:-5]+"_meta.yml"
            else:
                raise RuntimeError("Unknown output suffix: "+ofname)
            #
            print("Writing file: "+dfile)
            with open(dfile, 'w') as OUTPUT:
                yaml.dump(metadata, OUTPUT)
        #
        # Launch the task that is named by
        #
        tasks[args.subparser_name].run(
                    {},
                    cargs
                    )
