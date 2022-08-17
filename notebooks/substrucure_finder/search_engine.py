from pathlib import Path
from typing import Generator
from concurrent.futures import ThreadPoolExecutor, as_completed

from substrucure_finder import utils
from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import consts
from substrucure_finder.splitter_tree import SplitterTree
from substrucure_finder.bucket_search_engine import BucketSearchEngine


class SearchEngine:
    def __init__(self, data_path: str) -> None:
        self.data_path = Path(data_path)
        self.splitter_tree = SplitterTree.load(utils.splitter_tree_path(self.data_path))

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            yield from BucketSearchEngine.search_in_drive(fingerprint, self.data_path, bucket)


class ScopeExecutor(ThreadPoolExecutor):
    def __del__(self):
        print('del executor')
        self.shutdown(wait=False, cancel_futures=True)


class ThreadPoolSearchEngine(SearchEngine):
    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        executor = ScopeExecutor()
        buckets = self.splitter_tree.get_buckets(fingerprint)

        futures = [executor.submit(BucketSearchEngine.search_in_drive, fingerprint, self.data_path, bucket)
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
