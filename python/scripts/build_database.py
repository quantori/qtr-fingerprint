import argparse

import build_bucket_search_engines
import build_splitter_tree


def create_parser() -> argparse.ArgumentParser:
    splitter_tree_parser = build_splitter_tree.create_parser()
    bucket_search_engines_parser = build_bucket_search_engines.create_parser()
    parser = argparse.ArgumentParser(parents=[splitter_tree_parser, bucket_search_engines_parser],
                                     conflict_handler='resolve')
    return parser


def check_arguments(arguments: argparse.Namespace):
    build_splitter_tree.check_arguments(arguments)
    build_bucket_search_engines.check_arguments(arguments)


def build_database(arguments: argparse.Namespace):
    build_splitter_tree.build_splitter_tree(arguments)
    build_bucket_search_engines.build_bucket_search_engines(arguments)


if __name__ == '__main__':
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    check_arguments(args)
    build_database(args)
