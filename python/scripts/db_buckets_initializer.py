import argparse
import sys

sys.path.append('../packages/substructure_finder')
sys.path.append('../packages/')

from packages.substructure_finder import BucketsInitializerx`


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description='Init buckets from raw data base')
    parser.add_argument("--bt_dims", help="Ball tree dimensions number", type=int, required=True)
    parser.add_argument("--raw_db_name", help="Name of source raw data base", type=str, required=True)
    parser.add_argument("--db_name", help="Name of destination datat base", type=str, required=True)
    parser.add_argument("--data_paths", help="Paths to buckets data", nargs='+', type=str, required=True)
    parser.add_argument("--other_data_path", help="Path to other data", type=str, required=True)
    arguments = parser.parse_args(argv)
    return arguments


if __name__ == '__main__':
    args = parse_arguments()
