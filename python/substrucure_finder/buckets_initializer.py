from pathlib import Path
import pandas as pd
from pandarallel import pandarallel

from substrucure_finder.bucket_search_engine import BucketSearchEngine
from substrucure_finder import utils
from substrucure_finder.loader import Loader
from substrucure_finder.dumper import Dumper


class BucketsInitializer:
    def __init__(self, load_path: str, dump_path: str):
        self.load_path = Path(load_path)
        self.dump_path = Path(dump_path)

    def init_buckets(self) -> None:
        with Loader(utils.splitter_tree_path(self.load_path)) as loader:
            splitter_tree = loader.splitter_tree()
        pandarallel.initialize()
        buckets = pd.Series(list(splitter_tree.all_buckets))
        buckets.parallel_apply(self.__init_bucket)

    def __init_bucket(self, bucket: int) -> None:
        with Loader(utils.raw_bucket_path(self.load_path, bucket)) as loader:
            raw_bucket = loader.raw_bucket()
        with Loader(utils.columns_path(self.load_path, bucket)) as loader:
            columns = loader.columns()
        bucket_search_engine = BucketSearchEngine(raw_bucket, columns)
        with Dumper(utils.bucket_path(self.dump_path, bucket)) as dumper:
            dumper.bucket_search_engine(bucket_search_engine)
