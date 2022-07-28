from __future__ import annotations
from typing import List, Generator

from substrucure_finder import consts
from substrucure_finder.consts import Fingerprint


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

    def __walk_tree(self, fingerprint: Fingerprint, node: SplitterTree.Node) -> Generator[int, None, None]:
        if node.is_leaf:
            yield node.id
        else:
            if fingerprint[node.split_bit] == 0:
                yield from self.__walk_tree(fingerprint, self.nodes[node.left])
            yield from self.__walk_tree(fingerprint, self.nodes[node.right])

    def get_buckets(self, fingerprint: Fingerprint) -> Generator[int, None, None]:
        return self.__walk_tree(fingerprint, self.root)
