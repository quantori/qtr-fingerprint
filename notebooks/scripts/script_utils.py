import sys
import os


def run_command(command: str):
    print("Run command: ", command)
    code = os.system(command)
    if code != 0:
        print(f'Command:\n\t{command}\nfiled with exit code {code}')
        sys.exit(code)
