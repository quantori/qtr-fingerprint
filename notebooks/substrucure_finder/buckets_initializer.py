from pathlib import Path
import pandas as pd
from pandarallel import pandarallel

from substrucure_finder.bucket_search_engine import BucketSearchEngine
from substrucure_finder.splitter_tree import SplitterTree
from substrucure_finder.molecule import Molecule
from substrucure_finder import utils


class BucketsInitializer:
    def __init__(self, raw_db_path: Path, db_path: Path, columns_count: int):
        self.raw_db_path = raw_db_path
        self.db_path = db_path
        self.columns_count = columns_count

    def init_buckets(self) -> None:
        splitter_tree = SplitterTree.load(utils.splitter_tree_path(self.raw_db_path))
        pandarallel.initialize()
        buckets = pd.Series(list(sorted(splitter_tree.all_buckets)))
        # TODO parallel apply. pickle does not support multi-threads???
        buckets.apply(self.__init_bucket)

    def __init_bucket(self, bucket: int) -> None:
        print(f"Start init {bucket}")
        molecules_df = Molecule.load_list_from_dir(utils.raw_bucket_path(self.raw_db_path, bucket))
        columns = utils.load_columns_from_file(utils.columns_path(self.raw_db_path, bucket))[:self.columns_count]
        assert len(columns) == self.columns_count
        BucketSearchEngine(molecules_df, columns, self.db_path, bucket)
        print(f"Finish init {bucket}")
