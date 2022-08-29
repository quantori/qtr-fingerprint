import argparse
import shutil
import sys
from pathlib import Path

packages_dir = Path(__file__).absolute().parent.parent / 'packages'
assert packages_dir.is_dir()
sys.path.append(str(packages_dir))

# noinspection PyUnresolvedReferences
from substructure_finder import DbFilesystem
# noinspection PyUnresolvedReferences
from substructure_finder.buckets_initializer import BucketsInitializer
# noinspection PyUnresolvedReferences
from fp_utils import CatchTime


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_db_name', type=str, required=True)
    parser.add_argument('--db_name', type=str, required=True)
    parser.add_argument('--ball_tree_dims', type=int, required=True)
    parser.add_argument('--store_dir_paths', type=str, nargs='+', required=True)
    parser.add_argument('--other_data_path', type=str, required=True)
    return parser


def check_arguments(arguments: argparse.Namespace):
    assert all(map(lambda x: Path(x).is_dir(), arguments.store_dir_paths)), "Every store_dir_path must be a directory"
    assert Path(arguments.other_data_path).is_dir(), "other_data_path must be a directory"
    assert arguments.ball_tree_dims > 1, "ball_tree_dims must be greater than 1"


def build_bucket_search_engines(arguments: argparse.Namespace):
    db_fs = DbFilesystem(arguments.store_dir_paths, arguments.other_data_path)

    print(f'Copy splitter tree from {arguments.raw_db_name} to {arguments.db_name}')
    db_tree_path = db_fs.tree_path(arguments.db_name)
    db_tree_path.parent.mkdir(parents=True, exist_ok=True)
    with db_fs.tree_path(arguments.raw_db_name).open('rb') as fp_in:
        with db_tree_path.open('wb') as fp_out:
            shutil.copyfileobj(fp_in, fp_out)

    initializer = BucketsInitializer(db_fs, arguments.raw_db_name, arguments.db_name,
                                     columns_count=arguments.ball_tree_dims)

    catch_time = CatchTime("buckets initialization")
    with catch_time:
        initializer.init_buckets()
    print(catch_time)


if __name__ == '__main__':
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    check_arguments(args)
    build_bucket_search_engines(args)
