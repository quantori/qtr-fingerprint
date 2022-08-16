import asyncio
from pathlib import Path
from typing import Generator, Iterable, AsyncGenerator, Any
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

    def __del__(self):
        if hasattr(self, 'executor'):
            self.executor.shutdown()

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            yield from self.__search_in_bucket(fingerprint, bucket)

    def async_search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        buckets = self.splitter_tree.get_buckets(fingerprint)
        tasks = [asyncio.create_task(self.__async_search_in_bucket(fingerprint, bucket)) for bucket in buckets]
        for task in tasks:
            yield from await task

    async def __async_search_in_bucket(self, fingerprint: Fingerprint, bucket: int) -> Iterable[str]:
        bucket_search_engine = BucketSearchEngine.load(utils.bucket_path(self.data_path, bucket))
        return bucket_search_engine.search(fingerprint, self.data_path, bucket)

    def thread_pool_search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        if not hasattr(self, 'executor'):
            print('init executor')
            self.executor = ThreadPoolExecutor()
        buckets = self.splitter_tree.get_buckets(fingerprint)
        futures = [self.executor.submit(SearchEngine.__search_in_bucket, self, fingerprint, bucket)
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()

    def __search_in_bucket(self, fingerprint: Fingerprint, bucket: int) -> Iterable[str]:
        bucket_search_engine = BucketSearchEngine.load(utils.bucket_path(self.data_path, bucket))
        return bucket_search_engine.search(fingerprint, self.data_path, bucket)
