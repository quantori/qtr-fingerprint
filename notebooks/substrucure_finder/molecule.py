from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from typing import BinaryIO, List

from substrucure_finder import consts
from substrucure_finder.fingerprint import ByteFingerprint


@dataclass
class Molecule:
    byte_fingerprint: ByteFingerprint
    smiles: str

    @classmethod
    def load_byte_fingerprint_from_stream(cls, stream: BinaryIO) -> ByteFingerprint:
        fingerprint_bytes = stream.read(consts.fingerprint_size_in_bytes)
        assert len(fingerprint_bytes) == consts.fingerprint_size_in_bytes
        return ByteFingerprint(np.fromiter(fingerprint_bytes, dtype=np.uint8))

    @classmethod
    def load_smiles_from_stream(cls, stream: BinaryIO) -> str:
        return stream.readline().decode(encoding='utf-8').strip()

    @classmethod
    def load_molecule_from_stream(cls, stream: BinaryIO) -> Molecule:
        fp = cls.load_byte_fingerprint_from_stream(stream)
        smiles = cls.load_smiles_from_stream(stream)
        return Molecule(fp, smiles)

    @classmethod
    def load_molecules_from_rb_stream(cls, stream: BinaryIO) -> List[Molecule]:
        molecules_number = int.from_bytes(stream.read(8), byteorder=consts.byteorder, signed=False)
        result = []
        for _ in range(molecules_number):
            result.append(cls.load_molecule_from_stream(stream))
        return result
