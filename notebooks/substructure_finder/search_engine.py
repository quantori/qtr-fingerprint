from typing import Generator
from concurrent.futures import ThreadPoolExecutor, as_completed

from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.fingerprint import BitFingerprint
from substructure_finder import consts
from substructure_finder.splitter_tree import SplitterTree
from substructure_finder.bucket_search_engine import BucketSearchEngine


class SearchEngine:
    def __init__(self, db_filesystem: DbFilesystem, db_name: str) -> None:
        self.splitter_tree = SplitterTree.load(db_filesystem.tree_path(db_name))
        self.pickle_paths = db_filesystem.pickle_paths(db_name)

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            yield from BucketSearchEngine.search_in_file(fingerprint, self.pickle_paths[bucket])


class ScopeExecutor(ThreadPoolExecutor):
    def __del__(self):
        self.shutdown(wait=False, cancel_futures=True)


class ThreadPoolSearchEngine(SearchEngine):

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        executor = ScopeExecutor()
        buckets = self.splitter_tree.get_buckets(fingerprint)
        futures = [executor.submit(BucketSearchEngine.search_in_file, fingerprint, self.pickle_paths[bucket])
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
