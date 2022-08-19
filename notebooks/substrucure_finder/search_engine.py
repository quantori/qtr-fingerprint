import os
from pathlib import Path
from typing import Generator, List
from concurrent.futures import ThreadPoolExecutor, as_completed

from substrucure_finder import utils
from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import consts
from substrucure_finder.splitter_tree import SplitterTree
from substrucure_finder.bucket_search_engine import BucketSearchEngine


class SearchEngine:
    def __init__(self, data_paths: List[str], other_data_path: Path, db_name: str) -> None:
        self.data_paths = [Path(data_path) / db_name for data_path in data_paths]
        assert all(path.is_dir() for path in self.data_paths)
        self.other_data_path = Path(other_data_path) / db_name
        assert self.other_data_path.is_dir()
        self.splitter_tree = SplitterTree.load(utils.splitter_tree_path(self.other_data_path))
        self.buckets_paths = dict()
        all_buckets_names = set(self.splitter_tree.all_buckets)
        for buckets_data_path in self.data_paths:
            for bucket_name in map(str, os.listdir(buckets_data_path)):
                bucket_id = int(bucket_name)
                assert bucket_id in all_buckets_names, f"{bucket_name} is not a bucket name"
                bucket_path = buckets_data_path / bucket_name
                self.buckets_paths[bucket_id] = bucket_path
        assert len(all_buckets_names) >= len(self.buckets_paths), "No all buckets was found in data"
        assert len(all_buckets_names) <= len(self.buckets_paths), "Extra buckets was found in data"

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            yield from BucketSearchEngine.search_in_drive(fingerprint, self.buckets_paths[bucket], bucket)


class ScopeExecutor(ThreadPoolExecutor):
    def __del__(self):
        print('del executor')
        self.shutdown(wait=False, cancel_futures=True)


class ThreadPoolSearchEngine(SearchEngine):

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        executor = ScopeExecutor()
        buckets = self.splitter_tree.get_buckets(fingerprint)
        futures = [executor.submit(BucketSearchEngine.search_in_drive, fingerprint, self.buckets_paths[bucket], bucket)
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
