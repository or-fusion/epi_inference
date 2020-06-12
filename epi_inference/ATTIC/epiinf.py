#
# Main script for epi_inference
#
import argparse
from . import collect
from . import inference
from . import reconstruct

def main():
    Parser = argparse.ArgumentParser(description='inference models')
    Parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Verbosity flag')
    subparsers = Parser.add_subparsers(title='subcommands', help="Help", dest="subparser_name")
    #
    # epiinf collect
    #
    parser1 = subparsers.add_parser('collect', help='Collect data for inference models')
    parser1.add_argument('config_file',
                    help='YAML configuration file')
    parser1.set_defaults(subparser_func=collect.run)
    #
    # epiinf inference
    #
    parser2 = subparsers.add_parser('inference', help='Optimize an inference model')
    parser2.add_argument('config_file',
                    help='YAML configuration file')
    parser2.set_defaults(subparser_func=inference.run)
    #
    # epiinf reconstruct
    #
    parser3 = subparsers.add_parser('reconstruct', help='Reconstruct disease propigation given case data')
    parser3.add_argument('config_file',
                    help='YAML configuration file')
    parser3.set_defaults(subparser_func=reconstruct.run)
    #
    # Process arguments
    #
    args = Parser.parse_args()
    if args.subparser_name is None:
        Parser.print_help()
        return
    return args.subparser_func(args)


if __name__ == "__main__":
    main()

    

