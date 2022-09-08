from __future__ import annotations

from typing import List, Iterable

from sklearn.neighbors import BallTree
from pathlib import Path
import joblib

from substructure_finder import consts
from substructure_finder.fingerprint import BitFingerprint
from substructure_finder import utils
from substructure_finder.molecule import Molecule


class BucketSearchEngine:
    def __init__(self, molecules: List[Molecule], columns: List[int]) -> None:
        full_bit_fingerprints = [utils.byte_fingerprint_to_bits(mol.byte_fingerprint, consts.fingerprint_size_in_bits)
                                 for mol in molecules]
        sub_bit_fingerprints = [utils.take_columns_from_bit_fingerprint(columns, fp) for fp in full_bit_fingerprints]

        self.ball_tree = BallTree(sub_bit_fingerprints, leaf_size=2, metric='russellrao')
        self.columns = columns
        self.molecules = molecules

    def search(self, bit_fingerprint: BitFingerprint) -> Iterable[str]:
        assert bit_fingerprint.dtype == bool
        sub_fingerprint = utils.take_columns_from_bit_fingerprint(self.columns, bit_fingerprint)
        radius = utils.russelrao_radius(sub_fingerprint)
        ball_tree_answers = self.ball_tree.query_radius([sub_fingerprint], radius)[0]
        if len(ball_tree_answers) == 0:
            return list()

        bytes_fingerprint = utils.bit_fingerprint_to_bytes(bit_fingerprint)
        filtered_answers = [i for i in sorted(ball_tree_answers) if
                            utils.is_sub_fingerprint(bytes_fingerprint, self.molecules[i].byte_fingerprint)]
        if len(filtered_answers) == 0:
            return list()

        smiles_answers = [self.molecules[i].smiles for i in filtered_answers]
        return smiles_answers

    @classmethod
    def search_in_file(cls, fingerprint: BitFingerprint, file: Path) -> Iterable[str]:
        bucket_search_engine = cls.load(file)
        return bucket_search_engine.search(fingerprint)

    def dump(self, file_path: Path) -> None:
        with file_path.open('wb') as f:
            joblib.dump(self, f, protocol=-1)

    @classmethod
    def load(cls, file: Path) -> BucketSearchEngine:
        with file.open('rb') as stream:
            bucket_search_engine = joblib.load(stream)
            assert isinstance(bucket_search_engine, cls)
        return bucket_search_engine
