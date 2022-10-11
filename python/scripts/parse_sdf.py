import argparse
import os
import gzip
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
import requests
import pandas as pd
from pathlib import Path
from bs4 import BeautifulSoup

sys.path.append(str(Path(__file__).absolute().parent.parent))

from scripts import build_and_run

url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
processed_file = Path(__file__).absolute().parent / '../../data/pubchem/processed.txt'
sdf_dir_path = "not initialized path"
dest_dir_path = "not initialized path"
parse_mode = "not initialized mode"

target = "parse_sdf"


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf_dir_path", type=str, required=True)
    parser.add_argument("--dest_dir_path", type=str, required=True)
    parser.add_argument("--parse_mode", type=str, required=True)
    return parser


def get_urls(url, ext):
    page = requests.get(url)
    soup = BeautifulSoup(page.content, features='lxml')
    all_links = [link.get("href") for link in soup("a")]
    return pd.Series(filter(lambda x: x.endswith(ext), all_links))


def tag_processed(file_name):
    with open(processed_file, 'a') as f:
        f.write(str(file_name) + '\n')
    print(f'processed: {file_name}', file=sys.stderr)


def check_processed(file_name, processed):
    for processed_name in processed:
        if file_name in processed_name:
            return True
    return False


def remove_unprocessed():
    with processed_file.open('r') as f:
        processed = f.read().split()
    for file in os.listdir(sdf_dir_path):
        if check_processed(file, processed):
            continue
        if (sdf_dir_path / file).is_dir():
            os.rmdir(sdf_dir_path / file)
        else:
            os.remove(sdf_dir_path / file)
        print('remove:', sdf_dir_path / file, file=sys.stderr)


def download_zip(file_name, where_to):
    if file_name in processed_before:
        print(f"skipped downloading of {file_name}, because it's already processed", file=sys.stderr)
        return
    print(f"start download {file_name}", file=sys.stderr)
    with requests.get(url + file_name, stream=True) as r:
        with open(where_to / file_name, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    out_file_path = where_to / file_name[:-3]
    with gzip.open(where_to / file_name, 'rb') as f_in:
        with open(out_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(where_to / file_name)
    return out_file_path


def process(file_name):
    where_to = sdf_dir_path / str(file_name).split('.')[0]
    where_to.mkdir(exist_ok=True)
    sdf_file = download_zip(file_name, where_to)
    build_and_run.run(target, [f"--source_dir_path={sdf_file.parent}",
                               f"--dest_dir_path={dest_dir_path}",
                               f"--parse_mode={parse_mode}"])
    os.remove(sdf_file)
    os.rmdir(where_to)
    tag_processed(file_name)


if __name__ == "__main__":
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    print(args, file=sys.stderr)
    build_and_run.build(target)

    sdf_dir_path = Path(args.sdf_dir_path)
    dest_dir_path = Path(args.dest_dir_path)
    parse_mode = args.parse_mode

    remove_unprocessed()
    dest_dir_path.mkdir(parents=True, exist_ok=True)
    sdf_dir_path.mkdir(parents=True, exist_ok=True)
    with processed_file.open('a'):  # ensure that file exists
        pass
    with open(processed_file, 'r') as f:
        processed_before = set(f.read().split())
    sdfs = get_urls(url, ".sdf.gz")
    file_names = pd.Series(list(sorted(sdfs)))[:2]

    ThreadPoolExecutor().map(process, file_names)
