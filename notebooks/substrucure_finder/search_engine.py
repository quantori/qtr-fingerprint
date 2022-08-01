from pathlib import Path
from typing import Generator

from substrucure_finder.loader import Loader
from substrucure_finder import utils
from substrucure_finder.consts import Fingerprint
from substrucure_finder import consts


class SearchEngine:
    def __init__(self, data_path: str) -> None:
        self.data_path = Path(data_path)
        with Loader(utils.splitter_tree_path(self.data_path)) as loader:
            # TODO check splitter tree correctness
            self.splitter_tree = loader.splitter_tree()

    def search(self, fingerprint: Fingerprint) -> Generator[str, None, None]:
        assert len(fingerprint) == consts.fingerprint_size
        # TODO get_buckets(fingerprint)
        for bucket in self.splitter_tree.all_buckets:
            print(f'Look into bucket: {bucket}')
            with Loader(utils.bucket_path(self.data_path, bucket)) as loader:
                bucket_search_engine = loader.bucket_search_engine()
            assert set(bucket_search_engine.columns) == set(range(consts.fingerprint_size))
            yield from bucket_search_engine.search(fingerprint)
