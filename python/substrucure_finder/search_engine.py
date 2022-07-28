from pathlib import Path
from typing import Generator

from substrucure_finder.loader import Loader
from substrucure_finder import utils
from substrucure_finder.consts import Fingerprint


class SearchEngine:
    def __init__(self, data_path: str) -> None:
        self.data_path = Path(data_path)
        with Loader(utils.splitter_tree_path(self.data_path)) as loader:
            self.splitter_tree = loader.splitter_tree()

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            with Loader(self.data_path / utils.bucket_path(self.data_path, bucket)) as loader:
                bucket_search_engine = loader.bucket_search_engine()
            yield from bucket_search_engine.search(fingerprint)
