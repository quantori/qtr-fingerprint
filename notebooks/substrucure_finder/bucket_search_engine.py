from __future__ import annotations

from typing import List, Iterable

import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import joblib

from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import utils


class BucketSearchEngine:
    def __init__(self, molecules: pd.DataFrame, columns: List[int]) -> None:
        self.columns = columns
        self.molecules = molecules
        sub_fingerprints = [utils.take_columns_from_fingerprint(self.columns, fp) for fp in
                            map(lambda i: molecules.iloc[i].values, range(len(molecules)))]
        self.ball_tree = BallTree(sub_fingerprints, leaf_size=2, metric='russellrao')

    def search(self, fingerprint: Fingerprint, data_path: Path, bucket_id: int) -> Iterable[str]:
        assert fingerprint.dtype == bool
        sub_fingerprint = utils.take_columns_from_fingerprint(self.columns, fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        if len(ball_tree_answers) == 0:
            return list()

        filtered_answers = [i for i in ball_tree_answers if
                            utils.is_sub_fingerprint(fingerprint, self.molecules.iloc[i].values)]
        answers = [self.molecules.iloc[i].name for i in filtered_answers]
        return answers

    @classmethod
    def load(cls, file_path: Path) -> BucketSearchEngine:
        assert file_path.is_file(), f"Path to load bucket search engine from must be a file, got {file_path}"
        with file_path.open('rb') as stream:
            obj = joblib.load(stream)
        assert isinstance(obj, BucketSearchEngine)
        return obj

    def dump(self, file_path: Path) -> None:
        with file_path.open('wb') as f:
            joblib.dump(self, f, protocol=-1)
