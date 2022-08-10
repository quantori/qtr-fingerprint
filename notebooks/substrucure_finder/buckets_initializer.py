from pathlib import Path
import pandas as pd
from pandarallel import pandarallel

from substrucure_finder.bucket_search_engine import BucketSearchEngine
from substrucure_finder.splitter_tree import SplitterTree
from substrucure_finder.molecule import Molecule
from substrucure_finder import utils
from substrucure_finder import consts


class BucketsInitializer:
    def __init__(self, dir_path: Path, columns_count: int):
        self.dir_path = dir_path
        self.columns_count = columns_count

    def init_buckets(self) -> None:
        splitter_tree = SplitterTree.load(utils.splitter_tree_path(self.dir_path))
        pandarallel.initialize()
        buckets = pd.Series(list(sorted(splitter_tree.all_buckets)))
        # TODO parallel apply. pickle does not support multi-threads???
        buckets.apply(self.__init_bucket)

    def __init_bucket(self, bucket: int) -> None:
        print(f"Start init {bucket}")
        raw_bucket = Molecule.load_list_from_dir(utils.raw_bucket_path(self.dir_path, bucket))
        assert all(map(lambda x: len(x.fingerprint) == consts.fingerprint_size, raw_bucket))
        columns = utils.load_columns_from_file(utils.columns_path(self.dir_path, bucket))[:self.columns_count]
        assert len(columns) == self.columns_count
        bucket_search_engine = BucketSearchEngine(raw_bucket, columns)
        bucket_search_engine.dump(utils.bucket_path(self.dir_path, bucket))
        # smiles_list = list(map(lambda x: x.smiles, raw_bucket))
        # fingerprints_list = list(map(lambda x: x.fingerprint, raw_bucket))
        # Molecule.dump_smiles_list(smiles_list, self.dir_path / 'buckets' / (str(bucket) + '.smiles'))
        # Molecule.dump_fingerprints_list(fingerprints_list, self.dir_path / 'buckets' / (str(bucket) + '.fp'))
        print(f"Finish init {bucket}")
