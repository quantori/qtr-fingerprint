from typing import Mapping, List
import pickle
import numpy as np
from pathlib import Path

from substrucure_finder.splitter_tree import SplitterTree
from substrucure_finder import consts
from substrucure_finder.bucket_search_engine import BucketSearchEngine
from substrucure_finder.consts import Fingerprint
from substrucure_finder import utils


class Loader:
    def __init__(self, file: Path):
        assert (file.exists(), str(file))
        self.file = file

    def __enter__(self):
        self.stream = self.file.open('rb')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stream.close()

    def __splitter_tree_node(self) -> SplitterTree.Node:
        bytes_array = self.stream.read(32)
        assert len(bytes_array) == 32
        node_id = int.from_bytes(bytes_array[0:8], byteorder=consts.byteorder, signed=False)
        split_bit = int.from_bytes(bytes_array[8:16], byteorder=consts.byteorder, signed=False)
        left = int.from_bytes(bytes_array[16:24], byteorder=consts.byteorder, signed=False)
        right = int.from_bytes(bytes_array[24:32], byteorder=consts.byteorder, signed=False)
        return SplitterTree.Node(node_id, split_bit, left, right)

    def splitter_tree(self) -> SplitterTree:
        nodes_count = int.from_bytes(self.stream.read(8), byteorder=consts.byteorder, signed=False)
        nodes = [self.__splitter_tree_node() for _ in range(nodes_count)]
        nodes.sort(key=lambda node: node.id)
        return SplitterTree(nodes)

    def bucket_search_engine(self) -> BucketSearchEngine:
        obj = pickle.load(self.stream)
        assert isinstance(obj, BucketSearchEngine)
        return obj

    def fingerprint(self) -> Fingerprint:
        fingerprint_bytes = self.stream.read(consts.fingerprint_size_in_bytes)
        assert len(fingerprint_bytes) == consts.fingerprint_size_in_bytes
        fingerprint_bin_str = ''.join(map(utils.byte_to_bits, fingerprint_bytes))
        assert len(fingerprint_bin_str) == consts.fingerprint_size
        # TODO check byteorder and correctness of transformation above
        fingerprint = Fingerprint(np.array(list(map(int, fingerprint_bin_str))))
        assert len(fingerprint) == consts.fingerprint_size
        return fingerprint

    def smiles(self) -> str:
        symbols = []
        while True:
            symbol = self.stream.read(1)
            if len(symbol) == 0 or symbol[0] == ord('\n'):
                break
            symbols.append(symbol[0])
        return ''.join(map(lambda x: chr(x), symbols))

    def raw_bucket(self) -> Mapping[str, Fingerprint]:
        molecules_number = int.from_bytes(self.stream.read(8), byteorder=consts.byteorder, signed=False)
        result = dict()
        for _ in range(molecules_number):
            fingerprint = self.fingerprint()
            smiles = self.smiles()
            assert len(fingerprint) == consts.fingerprint_size
            result[smiles] = fingerprint
        return result

    def columns(self) -> List[int]:
        return list(map(int, str(self.stream.read().decode('utf-8')).split()))
