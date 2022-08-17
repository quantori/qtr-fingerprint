from __future__ import annotations

import pickle
from dataclasses import dataclass
from typing import List
from pathlib import Path
import numpy as np
import pandas as pd

from substrucure_finder import consts
from substrucure_finder.fingerprint import Fingerprint


@dataclass
class Molecule:
    smiles: str
    fingerprint: Fingerprint

    @classmethod
    def load_list_from_dir(cls, dir_path: Path) -> pd.DataFrame:
        file_path = dir_path / "0.csv"
        df = pd.read_csv(file_path, header=0, index_col=0, delimiter='~',
                         dtype=dict((str(i), bool) for i in range(consts.fingerprint_size_in_bits)))
        return df

    @classmethod
    def dump_fingerprints_list(cls, fingerprints: List[Fingerprint], file_path: Path) -> None:
        with file_path.open('wb') as stream:
            obj = np.concatenate(fingerprints, axis=0).astype(bool)
            pickle.dump(obj, stream, protocol=-1)

    @classmethod
    def load_fingerprints_list(cls, file_path: Path) -> List[Fingerprint]:
        with file_path.open('rb') as stream:
            obj = pickle.load(stream)
        return list(obj)

    @classmethod
    def dump_smiles_list(cls, smiles: List[str], file_path: Path) -> None:
        with file_path.open('wb') as stream:
            pickle.dump(smiles, stream, protocol=-1)

    @classmethod
    def load_smiles_list(cls, file_path: Path) -> List[str]:
        with file_path.open('rb') as stream:
            obj = pickle.load(stream)
        assert isinstance(obj, list)
        return obj
