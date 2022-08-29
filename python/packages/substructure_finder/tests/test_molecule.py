import unittest
import shutil
from pathlib import Path

import numpy as np

from substructure_finder import consts
from substructure_finder.fingerprint import ByteFingerprint
from substructure_finder.molecule import Molecule


class TestMolecule(unittest.TestCase):

    def setUp(self) -> None:
        self.tmp_dir = Path(__file__).absolute().parent / 'tmp'
        self.tmp_dir.mkdir()

        self.fp_byte_array = bytes([i % 256 for i in range(consts.fingerprint_size_in_bytes)])
        self.byte_fp = ByteFingerprint(np.fromiter(self.fp_byte_array, dtype=np.uint8))
        self.smiles = 'N-N-Cl'
        self.molecule = Molecule(self.byte_fp, self.smiles)

        self.fp_2_byte_array = bytes([123 for _ in range(consts.fingerprint_size_in_bytes)])
        self.byte_fp_2 = ByteFingerprint(np.fromiter(self.fp_2_byte_array, dtype=np.uint8))
        self.smiles_2 = 'C1=CC=C(C=C1)C=O'
        self.molecule_2 = Molecule(self.byte_fp_2, self.smiles_2)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_dir)

    def assertNpArraysEqual(self, arr_1: np.ndarray, arr_2: np.ndarray):
        self.assertTrue(np.all(arr_1 == arr_2), f'\nexpected: {arr_1}\n  actual: {arr_2}')

    def test_load_byte_fingerprint_from_file(self):
        file_path = self.tmp_dir / 'byte_fp.bin'
        with file_path.open('wb') as f:
            f.write(self.fp_byte_array)
        with file_path.open('rb') as f:
            actual_byte_fp = Molecule.load_byte_fingerprint_from_stream(f)
        self.assertNpArraysEqual(self.byte_fp, actual_byte_fp)

    def test_load_smiles_from_stream(self):
        file_path = self.tmp_dir / 'smiles.txt'
        with file_path.open('wb') as f:
            f.write(self.smiles.encode())
            f.write(b'\n')
        with file_path.open('rb') as f:
            actual_smiles = Molecule.load_smiles_from_stream(f)
        self.assertEqual(self.smiles, actual_smiles)

    def test_load_molecule_from_stream(self):
        file_path = self.tmp_dir / 'molecule.bin'
        with file_path.open('wb') as f:
            f.write(self.fp_byte_array)
            f.write(self.smiles.encode())
            f.write(b'\n')
        with file_path.open('rb') as f:
            molecule = Molecule.load_molecule_from_stream(f)
        self.assertNpArraysEqual(self.byte_fp, molecule.byte_fingerprint)
        self.assertEqual(self.smiles, molecule.smiles)

    def test_load_molecules_from_rb_stream(self):
        file_path = self.tmp_dir / 'data.rb'
        with file_path.open('wb') as f:
            f.write(int(2).to_bytes(length=8, byteorder=consts.byteorder, signed=False))
            f.write(self.fp_byte_array)
            f.write(self.smiles.encode())
            f.write(b'\n')
            f.write(self.fp_2_byte_array)
            f.write(self.smiles_2.encode())
            f.write(b'\n')
        with file_path.open('rb') as f:
            molecules = Molecule.load_molecules_from_rb_stream(f)

        self.assertEqual(2, len(molecules))
        self.assertEqual(self.smiles, molecules[0].smiles)
        self.assertEqual(self.smiles_2, molecules[1].smiles)
        self.assertNpArraysEqual(self.byte_fp, molecules[0].byte_fingerprint)
        self.assertNpArraysEqual(self.byte_fp_2, molecules[1].byte_fingerprint)
