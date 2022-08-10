from pathlib import Path
from typing import Generator
import sys

from substrucure_finder.loader import Loader
from substrucure_finder import utils
from substrucure_finder.consts import Fingerprint
from substrucure_finder import consts


class SearchEngine:
    def __init__(self, data_path: str) -> None:
        self.data_path = Path(data_path)
        self.splitter_tree = Loader(utils.splitter_tree_path(self.data_path)).splitter_tree()

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size
        for bucket in self.splitter_tree.get_buckets(fingerprint):
            # print(f'Look into bucket: {bucket}', file=sys.stderr)
            bucket_search_engine = Loader(utils.bucket_path(self.data_path, bucket)).bucket_search_engine()
            yield from bucket_search_engine.search(fingerprint)
