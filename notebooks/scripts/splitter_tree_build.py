from pathlib import Path
import argparse
from script_utils import run_command

project_dir = Path(__file__).absolute().parent.parent.parent
cpp_dir = project_dir / 'cpp'
build_dir = cpp_dir / "cmake-build-release"
target_name = 'splitter_tree'


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('--rb_dir_path', type=str, required=True)
    parser.add_argument('--store_dir_paths', type=str, nargs='+', required=True)
    parser.add_argument('--other_data_path', type=str, required=True)
    parser.add_argument('--cols_subset_path', type=str, required=True)
    parser.add_argument('--max_tree_depth', type=int, required=True)
    parser.add_argument('--stop_bucket_size', type=int, required=True)
    parser.add_argument('--subtree_parallel_depth', type=int, required=True)
    parser.add_argument('--raw_db_name', type=str, required=True)
    return parser


def check_arguments(arguments: argparse.Namespace) -> None:
    assert Path(arguments.rb_dir_path).is_dir(), "rb_dir_path must be a directory"
    assert all(map(lambda x: Path(x).is_dir(), arguments.store_dir_paths)), "Every store_dir_path must be a directory"
    assert Path(arguments.other_data_path).is_dir(), "other_data_path must be a directory"
    assert Path(arguments.cols_subset_path).is_file(), "cols_subset_path must be a file"
    assert arguments.max_tree_depth > 0, "max_tree_depth must be greater than 0"
    assert arguments.stop_bucket_size > 0, "stop_bucket_size must be grater than 0"
    assert arguments.subtree_parallel_depth > 0, "subtree_parallel_depth must be greater than 0"


def build_splitter_tree(arguments: argparse.Namespace):
    print('Init cmake:')
    init_cmake_command = f'cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=/usr/bin/ninja -G Ninja ' \
                         f'-S {cpp_dir} ' \
                         f'-B {build_dir}'
    run_command(init_cmake_command)

    print('Build project:')
    build_project_command = f'cmake --build {build_dir} --target {target_name} -j 6'
    run_command(build_project_command)

    print('Start splitter tree building:')
    build_splitter_tree_command = f'{build_dir / "bin" / target_name} ' \
                                  f'--rb_dir_path={arguments.rb_dir_path} ' \
                                  f'--store_dir_paths={",".join(arguments.store_dir_paths)} ' \
                                  f'--other_data_path={arguments.other_data_path} ' \
                                  f'--cols_subset_path={arguments.cols_subset_path} ' \
                                  f'--max_tree_depth={arguments.max_tree_depth} ' \
                                  f'--stop_bucket_size={arguments.stop_bucket_size} ' \
                                  f'--subtree_parallel_depth={arguments.subtree_parallel_depth} ' \
                                  f'--raw_db_name={arguments.raw_db_name}'
    run_command(build_splitter_tree_command)


if __name__ == "__main__":
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    check_arguments(args)
    build_splitter_tree(args)
