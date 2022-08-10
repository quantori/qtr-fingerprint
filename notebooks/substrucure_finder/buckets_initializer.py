from pathlib import Path
import pandas as pd
from pandarallel import pandarallel

from substrucure_finder.bucket_search_engine import BucketSearchEngine
from substrucure_finder import utils
from substrucure_finder.loader import Loader
from substrucure_finder.dumper import Dumper
from substrucure_finder import consts


class BucketsInitializer:
    def __init__(self, dir_path: Path, columns_count: int):
        self.dir_path = dir_path
        self.columns_count = columns_count

    def init_buckets(self) -> None:
        splitter_tree = Loader(utils.splitter_tree_path(self.dir_path)).splitter_tree()
        pandarallel.initialize()
        buckets = pd.Series(list(sorted(splitter_tree.all_buckets)))
        # TODO parallel apply. pickle does not support multi-threads???
        buckets.apply(self.__init_bucket)

    def __init_bucket(self, bucket: int) -> None:
        print(f"Start init {bucket}")
        raw_bucket = Loader(utils.raw_bucket_path(self.dir_path, bucket)).raw_bucket()
        assert all(map(lambda x: len(x) == consts.fingerprint_size, raw_bucket.values()))
        columns = Loader(utils.columns_path(self.dir_path, bucket)).columns()[:self.columns_count]
        assert len(columns) == self.columns_count
        bucket_search_engine = BucketSearchEngine(raw_bucket, columns)
        Dumper(utils.bucket_path(self.dir_path, bucket)).bucket_search_engine(bucket_search_engine)
        print(f"Finish init {bucket}")
