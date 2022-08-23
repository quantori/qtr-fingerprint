import sys
import unittest
import argparse
from pathlib import Path

import numpy as np

from substrucure_finder import utils
from substrucure_finder import consts
from substrucure_finder.fingerprint import BitFingerprint, ByteFingerprint


def assert_message(expected_val, actual_val) -> str:
    return f'\nexpected: {expected_val}\n  actual: {actual_val}'


class TestUtils(unittest.TestCase):

    def setUp(self) -> None:
        self.bit_fp_zeros = BitFingerprint(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=bool))
        self.bit_fp_ones = BitFingerprint(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=bool))
        self.bit_fp_1 = BitFingerprint(np.array([1, 0, 1, 0, 1, 0, 1, 0, 0], dtype=bool))
        self.bit_fp_2 = BitFingerprint(np.array([1, 1, 1, 0, 1, 1, 1, 0, 1], dtype=bool))
        self.bit_fp_3 = BitFingerprint(np.array([0, 1, 1, 0, 1, 0, 1, 0, 1], dtype=bool))

        self.byte_fp_zeros = ByteFingerprint(np.array([0, 0], dtype=np.uint8))
        self.byte_fp_ones = ByteFingerprint(np.array([255, 128], dtype=np.uint8))
        self.byte_fp_1 = ByteFingerprint(np.array([85, 0], dtype=np.uint8))
        self.byte_fp_2 = ByteFingerprint(np.array([119, 128], dtype=np.uint8))
        self.byte_fp_3 = ByteFingerprint(np.array([86, 128], dtype=np.uint8))

    def assertNpArraysEqual(self, arr_1: np.ndarray, arr_2: np.ndarray):
        self.assertTrue(np.all(arr_1 == arr_2), assert_message(arr_1, arr_2))

    def test_russelrao_radius(self):
        r_zeros = utils.russelrao_radius(self.bit_fp_zeros)
        self.assertAlmostEqual(1, r_zeros, delta=2 * consts.radius_eps)

        r_ones = utils.russelrao_radius(self.bit_fp_ones)
        self.assertAlmostEqual(0, r_ones, delta=2 * consts.radius_eps)

        r_1 = utils.russelrao_radius(self.bit_fp_1)
        self.assertAlmostEqual(5 / 9, r_1, delta=2 * consts.radius_eps)

        r_2 = utils.russelrao_radius(self.bit_fp_2)
        self.assertAlmostEqual(2 / 9, r_2, delta=2 * consts.radius_eps)

        r_3 = utils.russelrao_radius(self.bit_fp_3)
        self.assertAlmostEqual(4 / 9, r_3, delta=2 * consts.radius_eps)

    def test_take_columns_from_bit_fingerprint(self):
        sub_fp_1 = utils.take_columns_from_bit_fingerprint([0, 1, 3, 5, 7], self.bit_fp_1)
        self.assertTrue(np.all(sub_fp_1 == np.array([1, 0, 0, 0, 0], dtype=bool)))

        sub_fp_2 = utils.take_columns_from_bit_fingerprint([2, 3, 4], self.bit_fp_2)
        self.assertTrue(np.all(sub_fp_2 == np.array([1, 0, 1], dtype=bool)))

        sub_fp_3 = utils.take_columns_from_bit_fingerprint([7, 6, 1, 2], self.bit_fp_3)
        self.assertTrue(np.all(sub_fp_3 == np.array([0, 1, 1, 1], dtype=bool)))

    def test_is_sub_fingerprint_bits(self):
        self.assertTrue(utils.is_sub_fingerprint(self.bit_fp_zeros, self.bit_fp_zeros))
        self.assertTrue(utils.is_sub_fingerprint(self.bit_fp_1, self.bit_fp_1))

        self.assertTrue(utils.is_sub_fingerprint(self.bit_fp_1, self.bit_fp_2))
        self.assertFalse(utils.is_sub_fingerprint(self.bit_fp_1, self.bit_fp_3))
        self.assertFalse(utils.is_sub_fingerprint(self.bit_fp_2, self.bit_fp_1))
        self.assertFalse(utils.is_sub_fingerprint(self.bit_fp_2, self.bit_fp_3))
        self.assertFalse(utils.is_sub_fingerprint(self.bit_fp_3, self.bit_fp_1))
        self.assertTrue(utils.is_sub_fingerprint(self.bit_fp_3, self.bit_fp_2))

    def test_is_sub_fingerprint_bytes(self):
        self.assertTrue(utils.is_sub_fingerprint(self.byte_fp_zeros, self.byte_fp_zeros))
        self.assertTrue(utils.is_sub_fingerprint(self.byte_fp_1, self.byte_fp_1))

        self.assertTrue(utils.is_sub_fingerprint(self.byte_fp_1, self.byte_fp_2))
        self.assertFalse(utils.is_sub_fingerprint(self.byte_fp_1, self.byte_fp_3))
        self.assertFalse(utils.is_sub_fingerprint(self.byte_fp_2, self.byte_fp_1))
        self.assertFalse(utils.is_sub_fingerprint(self.byte_fp_2, self.byte_fp_3))
        self.assertFalse(utils.is_sub_fingerprint(self.byte_fp_3, self.byte_fp_1))
        self.assertTrue(utils.is_sub_fingerprint(self.byte_fp_3, self.byte_fp_2))

    def test_byte_to_bits(self):
        self.assertEqual('01010000', utils.byte_to_bits(10))
        self.assertEqual('11111111', utils.byte_to_bits(255))
        self.assertEqual('00000000', utils.byte_to_bits(0))
        self.assertEqual('00100010', utils.byte_to_bits(68))
        with self.assertRaises(AssertionError):
            utils.byte_to_bits(-1)
        with self.assertRaises(AssertionError):
            utils.byte_to_bits(256)

    def test_bits_to_byte(self):
        self.assertEqual(10, utils.bits_to_byte(np.array([0, 1, 0, 1, 0, 0, 0, 0])))
        self.assertEqual(255, utils.bits_to_byte(np.array([1, 1, 1, 1, 1, 1, 1, 1])))
        self.assertEqual(0, utils.bits_to_byte(np.array([0, 0, 0, 0, 0, 0, 0, 0])))
        self.assertEqual(68, utils.bits_to_byte(np.array([0, 0, 1, 0, 0, 0, 1, 0])))
        self.assertEqual(128, utils.bits_to_byte(np.array([1])))

    def test_bit_fingerprint_to_str(self):
        self.assertEqual('000000000', utils.bit_fingerprint_to_str(self.bit_fp_zeros))
        self.assertEqual('111111111', utils.bit_fingerprint_to_str(self.bit_fp_ones))
        self.assertEqual('101010100', utils.bit_fingerprint_to_str(self.bit_fp_1))
        self.assertEqual('111011101', utils.bit_fingerprint_to_str(self.bit_fp_2))
        self.assertEqual('011010101', utils.bit_fingerprint_to_str(self.bit_fp_3))

    def test_bit_fingerprint_to_bytes(self):
        self.assertNpArraysEqual(self.byte_fp_zeros, utils.bit_fingerprint_to_bytes(self.bit_fp_zeros))
        self.assertNpArraysEqual(self.byte_fp_ones, utils.bit_fingerprint_to_bytes(self.bit_fp_ones))
        self.assertNpArraysEqual(self.byte_fp_1, utils.bit_fingerprint_to_bytes(self.bit_fp_1))
        self.assertNpArraysEqual(self.byte_fp_2, utils.bit_fingerprint_to_bytes(self.bit_fp_2))
        self.assertNpArraysEqual(self.byte_fp_3, utils.bit_fingerprint_to_bytes(self.bit_fp_3))

    def test_byte_fingerprint_to_bits(self):
        self.assertNpArraysEqual(self.bit_fp_zeros, utils.byte_fingerprint_to_bits(self.byte_fp_zeros, 9))
        self.assertNpArraysEqual(self.bit_fp_ones, utils.byte_fingerprint_to_bits(self.byte_fp_ones, 9))
        self.assertNpArraysEqual(self.bit_fp_1, utils.byte_fingerprint_to_bits(self.byte_fp_1, 9))
        self.assertNpArraysEqual(self.bit_fp_2, utils.byte_fingerprint_to_bits(self.byte_fp_2, 9))
        self.assertNpArraysEqual(self.bit_fp_3, utils.byte_fingerprint_to_bits(self.byte_fp_3, 9))


if __name__ == '__main__':
    unittest.main()
