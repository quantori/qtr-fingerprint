from __future__ import annotations
from typing import List, Generator, BinaryIO
from pathlib import Path

from substructure_finder import consts
from substructure_finder.fingerprint import BitFingerprint


class SplitterTree:
    class Node:
        def __init__(self, node_id: int, split_bit: int, left: int, right: int) -> None:
            self.id = node_id
            self.split_bit = split_bit
            self.left = left
            self.right = right

        @property
        def is_leaf(self) -> bool:
            return self.split_bit == consts.uint64_minus_one

        @classmethod
        def load_from_stream(cls, stream: BinaryIO) -> SplitterTree.Node:
            bytes_array = stream.read(32)
            assert len(bytes_array) == 32
            node_id = int.from_bytes(bytes_array[0:8], byteorder=consts.byteorder, signed=False)
            split_bit = int.from_bytes(bytes_array[8:16], byteorder=consts.byteorder, signed=False)
            left = int.from_bytes(bytes_array[16:24], byteorder=consts.byteorder, signed=False)
            right = int.from_bytes(bytes_array[24:32], byteorder=consts.byteorder, signed=False)
            return SplitterTree.Node(node_id, split_bit, left, right)

    def __init__(self, nodes: List[SplitterTree.Node]) -> None:
        self.nodes = nodes
        assert all(self.nodes[i].id == i for i in range(len(nodes)))

    @property
    def root(self):
        return self.nodes[0]

    @property
    def nodes_number(self):
        return len(self.nodes)

    def __find_leafs(self, node: SplitterTree.Node) -> Generator[int, None, None]:
        if node.is_leaf:
            yield node.id
        else:
            yield from self.__find_leafs(self.nodes[node.left])
            yield from self.__find_leafs(self.nodes[node.right])

    @property
    def all_buckets(self) -> Generator[int, None, None]:
        return self.__find_leafs(self.root)

    def __walk_tree(self, fingerprint: BitFingerprint, node: SplitterTree.Node) -> Generator[int, None, None]:
        if node.is_leaf:
            yield node.id
        else:
            if fingerprint[node.split_bit] == 0:
                yield from self.__walk_tree(fingerprint, self.nodes[node.left])
            yield from self.__walk_tree(fingerprint, self.nodes[node.right])

    def get_buckets(self, fingerprint: BitFingerprint) -> Generator[int, None, None]:
        return self.__walk_tree(fingerprint, self.root)

    @classmethod
    def load(cls, file_path: Path) -> SplitterTree:
        assert file_path.is_file(), "Path to load splitter tree from must be a file"
        with file_path.open('rb') as stream:
            nodes_count = int.from_bytes(stream.read(8), byteorder=consts.byteorder, signed=False)
            nodes = [SplitterTree.Node.load_from_stream(stream) for _ in range(nodes_count)]
            nodes.sort(key=lambda node: node.id)
            return SplitterTree(nodes)
