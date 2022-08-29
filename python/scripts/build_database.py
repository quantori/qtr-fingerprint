import argparse

import bucket_search_engines_build
import splitter_tree_build


def create_parser() -> argparse.ArgumentParser:
    splitter_tree_parser = splitter_tree_build.create_parser()
    bucket_search_engines_parser = bucket_search_engines_build.create_parser()
    parser = argparse.ArgumentParser(parents=[splitter_tree_parser, bucket_search_engines_parser],
                                     conflict_handler='resolve')
    return parser


def check_arguments(arguments: argparse.Namespace):
    splitter_tree_build.check_arguments(arguments)
    bucket_search_engines_build.check_arguments(arguments)


def build_database(arguments: argparse.Namespace):
    splitter_tree_build.build_splitter_tree(arguments)
    bucket_search_engines_build.build_bucket_search_engines(arguments)


if __name__ == '__main__':
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    check_arguments(args)
    build_database(args)
