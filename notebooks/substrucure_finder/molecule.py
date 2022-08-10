from __future__ import annotations

import pickle
import sys
from dataclasses import dataclass
from typing import BinaryIO, List
from pathlib import Path
import numpy as np

from substrucure_finder.buffered_stream import BufferedStream
from substrucure_finder import consts
from substrucure_finder.fingerprint import Fingerprint
from substrucure_finder.utils import byte_to_bits


@dataclass
class Molecule:
    smiles: str
    fingerprint: Fingerprint

    @classmethod
    def load_fingerprint_from_stream(cls, stream: BufferedStream) -> Fingerprint:
        fingerprint_bytes = stream.read(consts.fingerprint_size_in_bytes)
        assert len(fingerprint_bytes) == consts.fingerprint_size_in_bytes
        fingerprint_bin_str = ''.join(map(byte_to_bits, fingerprint_bytes))
        assert len(fingerprint_bin_str) == consts.fingerprint_size
        fingerprint = Fingerprint(np.fromiter(list(map(int, fingerprint_bin_str)), dtype=bool))
        assert len(fingerprint) == consts.fingerprint_size
        return fingerprint

    @classmethod
    def load_smiles_from_stream(cls, stream: BufferedStream) -> str:
        return stream.read_until('\n')

    @classmethod
    def load_list_from_stream(cls, stream: BufferedStream) -> List[Molecule]:
        molecules_number = int.from_bytes(stream.read(8), byteorder=consts.byteorder, signed=False)
        result = list()
        for _ in range(molecules_number):
            fingerprint = cls.load_fingerprint_from_stream(stream)
            smiles = cls.load_smiles_from_stream(stream)
            assert len(fingerprint) == consts.fingerprint_size
            result.append(Molecule(smiles, fingerprint))
        return result

    @classmethod
    def load_list_from_dir(cls, dir_path: Path) -> List[Molecule]:
        # print("start info loading", file=sys.stderr)
        stream = BufferedStream(dir_path / '0.rb')
        res = cls.load_list_from_stream(stream)
        # print("finish info loading", file=sys.stderr)
        return res
        # with (dir_path / '0.rb').open('rb') as stream:
        #     result += cls.load_list_from_stream(stream)
        # return result
        # assert dir_path.is_dir(), "Path to load bucket from must be a directory"
        # for file_name in os.listdir(dir_path):
        #     file_path = dir_path / str(file_name)
        #     assert file_path.is_file(), "Path to raw bucket file must be a file"
        #     with file_path.open('rb') as stream:
        #         result += cls.load_list_from_stream(stream)
        # return result

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
