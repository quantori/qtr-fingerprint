import os
import sys
from pathlib import Path
import argparse

project_dir = Path(__file__).absolute().parent.parent.parent
cpp_dir = project_dir / 'cpp'
build_dir = cpp_dir / "cmake-build-release"
target_name = 'ball_tree'


def run_command(command: str):
    print("Run command: ", command)
    code = os.system(command)
    if code != 0:
        print(f'Command:\n\t{command}\nfailed with exit code {code}')
        sys.exit(code)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('--rb_dir_path', type=str, required=True)
    parser.add_argument('--data_dir_paths', type=str, nargs='+', required=True)
    parser.add_argument('--other_data_path', type=str, required=True)
    parser.add_argument('--tree_depth', type=int, required=True)
    parser.add_argument('--subtree_parallel_depth', type=int, required=True)
    parser.add_argument('--db_name', type=str, required=True)
    return parser


def check_arguments(arguments: argparse.Namespace) -> None:
    assert Path(arguments.rb_dir_path).is_dir(), f"rb_dir_path must be a directory. Given path: {arguments.rb_dir_path}"
    assert all(map(lambda x: Path(x).is_dir(), arguments.data_dir_paths)), \
        f"Every data_dir_path must be a directory. " \
        f"Given paths: {list(map(lambda x: str(Path(x).absolute()), arguments.data_dir_paths))}"
    assert Path(arguments.other_data_path).is_dir(), "other_data_path must be a directory"
    assert arguments.tree_depth > 0, "tree_depth must be greater than 0"
    assert arguments.subtree_parallel_depth > 0, "subtree_parallel_depth must be greater than 0"


def build_ball_tree(arguments: argparse.Namespace):
    print('Init cmake:')
    init_cmake_command = f'cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=/usr/bin/ninja -G Ninja ' \
                         f'-S {cpp_dir} ' \
                         f'-B {build_dir}'
    run_command(init_cmake_command)

    print('Build project:')
    build_project_command = f'cmake --build {build_dir} --target {target_name} -j 6'
    run_command(build_project_command)

    print('Start ball tree building:')
    build_ball_tree_command = f'{build_dir / "bin" / target_name} ' \
                              f'--rb_dir_path={arguments.rb_dir_path} ' \
                              f'--data_dir_paths={",".join(arguments.data_dir_paths)} ' \
                              f'--other_data_path={arguments.other_data_path} ' \
                              f'--tree_depth={arguments.tree_depth} ' \
                              f'--subtree_parallel_depth={arguments.subtree_parallel_depth} ' \
                              f'--db_name={arguments.db_name}'
    run_command(build_ball_tree_command)


if __name__ == "__main__":
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    check_arguments(args)
    build_ball_tree(args)
