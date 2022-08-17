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
        self.molecules_cnt = len(molecules)

        self.dump(self, full_fingerprints, molecules.index.to_numpy(), utils.bucket_path(data_path, bucket_id))

    def search(self, fingerprint: Fingerprint, data_stream: BinaryIO, stream_start: int) -> Iterable[str]:
        assert fingerprint.dtype == bool
        sub_fingerprint = utils.take_columns_from_fingerprint(self.columns, fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        if len(ball_tree_answers) == 0:
            return list()

        bytes_fingerprint = utils.fingerprint_to_bytes(fingerprint)
        filtered_answers = [i for i in sorted(ball_tree_answers) if
                            self.check_answer(bytes_fingerprint, data_stream, stream_start, i)]
        if len(filtered_answers) == 0:
            return list()

        smiles_data_start = stream_start + self.molecules_cnt * consts.fingerprint_size_in_bytes
        smiles = self.load_smiles(data_stream, smiles_data_start)
        smiles_answers = [smiles[i] for i in filtered_answers]
        return smiles_answers

    @classmethod
    def check_answer(cls, byte_query: np.ndarray, fp_steam: BinaryIO, fp_data_start: int, ans_id: int) -> bool:
        fp_steam.seek(ans_id * consts.fingerprint_size_in_bytes + fp_data_start)
        ans_fp = np.fromiter(fp_steam.read(consts.fingerprint_size_in_bytes), dtype=np.uint8)
        return utils.is_sub_bytes_fingerprint(byte_query, ans_fp)

    @classmethod
    def search_in_drive(cls, fingerprint: Fingerprint, data_path: Path, bucket: int) -> Iterable[str]:
        with utils.bucket_path(data_path, bucket).open('rb') as stream:
            bucket_search_engine = joblib.load(stream)
            assert isinstance(bucket_search_engine, cls)
            data_start_pos = stream.tell()
            return bucket_search_engine.search(fingerprint, stream, data_start_pos)

    @classmethod
    def load_smiles(cls, stream: BinaryIO, smiles_data_start: int) -> List[str]:
        stream.seek(smiles_data_start)
        return stream.read().decode(encoding='utf-8').split(consts.smiles_split_symbol)

    @classmethod
    def dump(cls, search_engine: BucketSearchEngine, full_fingerprints: np.ndarray, smiles: Iterable[str],
             file_path: Path) -> None:

        bytes_fingerprints = np.apply_along_axis(utils.fingerprint_to_bytes, 1, full_fingerprints).astype(np.uint8)

        byte_str = b''.join(map(bytes, bytes_fingerprints))
        smiles_str = consts.smiles_split_symbol.join(smiles).encode(encoding='utf-8')

        with file_path.open('wb') as f:
            joblib.dump(search_engine, f, protocol=-1)
            f.write(byte_str)
            f.write(smiles_str)
