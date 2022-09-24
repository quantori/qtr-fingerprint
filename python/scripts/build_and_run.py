import os
import sys
from pathlib import Path
from typing import List

project_dir = Path(__file__).absolute().parent.parent.parent
cpp_dir = project_dir / 'cpp'
build_dir = cpp_dir / "cmake-build-release"


def run_command(command: str):
    print("Run command: ", command)
    code = os.system(command)
    if code != 0:
        print(f'Command:\n\t{command}\nfailed with exit code {code}')
        sys.exit(code)


def build_and_run(target: str, arguments: List[str]):
    print(f'Build and run: {target}')
    print('Init cmake:')
    init_cmake_command = f'cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_MAKE_PROGRAM=/usr/bin/ninja -G Ninja ' \
                         f'-S {cpp_dir} ' \
                         f'-B {build_dir}'
    run_command(init_cmake_command)

    print('Build project:')
    build_project_command = f'cmake --build {build_dir} --target {target} -j 6'
    run_command(build_project_command)

    print('Start program:')
    start_command = f'{build_dir / "bin" / target} ' + ' '.join(arguments)
    run_command(start_command)


if __name__ == "__main__":
    assert len(sys.argv) > 1
    build_and_run(sys.argv[1], sys.argv[2::])
