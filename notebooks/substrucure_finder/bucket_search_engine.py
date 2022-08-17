from __future__ import annotations

from typing import List, Iterable, BinaryIO

import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from pathlib import Path
import joblib

from substrucure_finder import consts
from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder import utils


class BucketSearchEngine:
    def __init__(self, molecules: pd.DataFrame, columns: List[int], data_path: Path, bucket_id: int) -> None:
        full_fingerprints = molecules.to_numpy(dtype=bool)
        sub_fingerprints = [utils.take_columns_from_fingerprint(columns, fp) for fp in full_fingerprints]

        self.ball_tree = BallTree(sub_fingerprints, leaf_size=2, metric='russellrao')
        self.columns = columns

        utils.bucket_path(data_path, bucket_id).mkdir()
        self.dump(utils.bucket_search_engine_path(data_path, bucket_id))
        self.dump_fingerprints(full_fingerprints, utils.full_fingerprints_path(data_path, bucket_id))
        self.dump_smiles(molecules.index.to_numpy(), utils.smiles_path(data_path, bucket_id))

    def search(self, fingerprint: Fingerprint, data_path: Path, bucket_id: int) -> Iterable[str]:
        assert fingerprint.dtype == bool
        sub_fingerprint = utils.take_columns_from_fingerprint(self.columns, fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        if len(ball_tree_answers) == 0:
            return list()

        bytes_fingerprint = utils.fingerprint_to_bytes(fingerprint)
        with utils.full_fingerprints_path(data_path, bucket_id).open('rb') as fp_stream:
            filtered_answers = [i for i in ball_tree_answers if self.check_answer(bytes_fingerprint, fp_stream, i)]
        if len(filtered_answers) == 0:
            return list()

        smiles = self.load_smiles(utils.smiles_path(data_path, bucket_id))
        smiles_answers = [smiles[i] for i in filtered_answers]
        return smiles_answers

    @classmethod
    def search_in_drive(cls, fingerprint: Fingerprint, data_path: Path, bucket: int) -> Iterable[str]:
        bucket_search_engine = cls.load(utils.bucket_search_engine_path(Path(data_path), bucket))
        return bucket_search_engine.search(fingerprint, data_path, bucket)

    @classmethod
    def load(cls, file_path: Path) -> BucketSearchEngine:
        assert file_path.is_file(), f"Path to load bucket search engine from must be a file, got {file_path}"
        with file_path.open('rb') as stream:
            obj = joblib.load(stream)
        assert isinstance(obj, BucketSearchEngine)
        return obj

    @classmethod
    def load_smiles(cls, file_path: Path) -> List[str]:
        with file_path.open('r') as f:
            return f.read().split(consts.smiles_split_symbol)

    def dump(self, file_path: Path) -> None:
        with file_path.open('wb') as f:
            joblib.dump(self, f, protocol=-1)

    @classmethod
    def check_answer(cls, byte_query: np.ndarray, fp_steam: BinaryIO, ans_id: int) -> bool:
        fp_steam.seek(ans_id * consts.fingerprint_size_in_bytes)
        ans_fp = np.fromiter(fp_steam.read(consts.fingerprint_size_in_bytes), dtype=np.uint8)
        return utils.is_sub_bytes_fingerprint(byte_query, ans_fp)

    @classmethod
    def dump_fingerprints(cls, fingerprints: np.ndarray, file_path: Path) -> None:
        bytes_fingerprints = np.apply_along_axis(utils.fingerprint_to_bytes, 1, fingerprints).astype(np.uint8)
        assert bytes_fingerprints.dtype == np.uint8
        assert all(map(lambda fp: len(fp) == consts.fingerprint_size_in_bytes, bytes_fingerprints))
        byte_str = b''.join(map(bytes, bytes_fingerprints))
        with file_path.open('wb') as f:
            f.write(byte_str)

    @classmethod
    def dump_smiles(cls, smiles: Iterable[str], file_path: Path) -> None:
        smiles_str = consts.smiles_split_symbol.join(smiles)
        with file_path.open('w') as f:
            f.write(smiles_str)
