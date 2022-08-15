from pathlib import Path
from typing import Generator, Iterable
import asyncio

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
            bucket_search_engine = BucketSearchEngine.load(utils.bucket_path(self.data_path, bucket))
            yield from bucket_search_engine.search(fingerprint, self.data_path, bucket)

    async def __search_in_bucket(self, fingerprint: Fingerprint, bucket: int) -> Iterable[str]:
        bucket_search_engine = BucketSearchEngine.load(utils.bucket_path(self.data_path, bucket))
        return bucket_search_engine.search(fingerprint, self.data_path, bucket)

    def async_search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size_in_bits

        buckets = self.splitter_tree.get_buckets(fingerprint)
        tasks = [asyncio.create_task(self.__search_in_bucket(fingerprint, bucket) for bucket in buckets)]

        for task in tasks:
            yield from await task
