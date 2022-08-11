from __future__ import annotations

import sys
from typing import List, Iterable

import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import joblib
import gzip
import time

from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import utils
from substrucure_finder import consts
from substrucure_finder.molecule import Molecule


class BucketSearchEngine:
    def __init__(self, molecules: pd.DataFrame, columns: List[int]) -> None:
        self.columns = columns
        self.molecules = molecules
        sub_fingerprints = [utils.take_columns_from_fingerprint(self.columns, fp) for fp in
                            map(lambda i: molecules.iloc[i].values, range(len(molecules)))]
        self.ball_tree = BallTree(sub_fingerprints, leaf_size=2, metric='russellrao')

    def search(self, fingerprint: Fingerprint, data_path: Path, bucket_id: int) -> Iterable[str]:
        # print('ball tree query', file=sys.stderr)
        assert fingerprint.dtype == bool
        sub_fingerprint = utils.take_columns_from_fingerprint(self.columns, fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        if len(ball_tree_answers) == 0:
            # print('no answers in ball tree', file=sys.stderr)
            return list()

        filtered_answers = [i for i in ball_tree_answers if
                            utils.is_sub_fingerprint(fingerprint, self.molecules.iloc[i].values)]
        answers = [self.molecules.iloc[i].name for i in filtered_answers]
        return answers
        # print('load more info', file=sys.stderr)
        # molecules = Molecule.load_list_from_dir(utils.raw_bucket_path(data_path, bucket_id))
        # print('load full fingerprints', file=sys.stderr)
        # fingerprints = Molecule.load_fingerprints_list(data_path / 'buckets' / (str(bucket_id) + '.fp'))
        # filtered_answers = [i for i in ball_tree_answers if
        #                     utils.is_sub_fingerprint(fingerprint, fingerprints[i])]
        # if len(filtered_answers) == 0:
        #     return list()
        # print('load smiles', file=sys.stderr)
        # smiles_list = Molecule.load_smiles_list(data_path / 'buckets' / (str(bucket_id) + '.smiles'))
        # answers = [smiles_list[i] for i in filtered_answers]
        # print('return answers', file=sys.stderr)
        # return answers

    @classmethod
    def load(cls, file_path: Path) -> BucketSearchEngine:
        assert file_path.is_file(), "Path to load bucket search engine from must be a file"
        with file_path.open('rb') as stream:
            obj = joblib.load(stream)
        assert isinstance(obj, BucketSearchEngine)
        return obj

        # compression loading
        # with gzip.open(file_path, 'rb') as stream:
        #     obj = joblib.load(stream)
        #     assert isinstance(obj, BucketSearchEngine)
        #     return obj

    def dump(self, file_path: Path) -> None:
        with file_path.open('wb') as f:
            joblib.dump(self, f, protocol=-1)

        # compression dumping:
        # print("start compression", file=sys.stderr)
        # p1 = time.time()
        # with gzip.open(file_path, 'wb', compresslevel=1) as f:
        #     joblib.dump(self, f, protocol=-1)
        # p2 = time.time()
        # print(f"finish compression {p2 - p1}", file=sys.stderr)
