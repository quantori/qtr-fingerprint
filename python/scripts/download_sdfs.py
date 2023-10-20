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

url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
dest_dir = "not initialized path"


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest_dir", type=str, required=True)
    return parser


def get_urls(url, ext):
    page = requests.get(url)
    soup = BeautifulSoup(page.content, features='lxml')
    all_links = [link.get("href") for link in soup("a")]
    return pd.Series(filter(lambda x: x.endswith(ext), all_links))


def download_zip(file_name, where_to):
    print(f"start download {file_name} to {where_to / file_name}", file=sys.stderr)
    with requests.get(url + file_name, stream=True) as r:
        with open(where_to / file_name, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def process(file_name):
    download_zip(file_name, dest_dir)


if __name__ == "__main__":
    arguments_parser = create_parser()
    args = arguments_parser.parse_args()
    print(args, file=sys.stderr)

    dest_dir = Path(args.dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    sdfs = get_urls(url, ".sdf.gz")
    file_names = pd.Series(list(sorted(sdfs)))

    ThreadPoolExecutor().map(process, file_names)
