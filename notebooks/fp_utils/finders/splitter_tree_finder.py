from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import List, Iterable, Optional, Tuple, Generator, Type
from sklearn.neighbors import BallTree
from abc import abstractmethod
from fp_utils.consts import PathType
from fp_utils.finders.drive_finder import DriveFinder
from fp_utils.finders.finder import Finder
from fp_utils.logger import Logger


class SplitterTreeFinder(DriveFinder):
    split_value = 0.5
    leaf_id = -1

    max_depth = 100

    @property
    @abstractmethod
    def bucket_size(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def inner_finder_class(self) -> Type[Finder]:
        raise NotImplementedError

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        return self.__tree.find(fingerprint, '')

    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        self.finder_path = Path(directory) / finder_name
        self.finder_path.mkdir(parents=True, exist_ok=True)
        self.__tree = self.__build_tree(df, list(), self.max_depth, '', *args, **kwargs)

    class Tree:
        def __init__(self, base: SplitterTreeFinder, col_id: int, left: Optional[SplitterTreeFinder.Tree],
                     right: Optional[SplitterTreeFinder.Tree], df: pd.DataFrame,
                     tree_path: str, *args, **kwargs) -> None:
            self.base = base
            self.col_id = col_id
            self.left = left
            self.right = right
            if self.is_leaf:
                self.inner_finder = self.base.inner_finder_class(df=df, directory=self.base.finder_path,
                                                                 finder_name=tree_path, *args, **kwargs)

        @property
        def is_leaf(self) -> bool:
            return self.col_id == self.base.leaf_id

        def find(self, fingerprint: pd.Series, node_path: str) -> Generator[str, None, None]:
            if self.is_leaf:
                yield from self.inner_finder.find_all(fingerprint)
            else:
                if fingerprint[self.col_id] == 0:
                    yield from self.left.find(fingerprint, node_path + '0')
                yield from self.right.find(fingerprint, node_path + '1')

    def __split_df(self, df: pd.DataFrame, used_cols: List[int]) -> Tuple[pd.Index, pd.Index, int]:
        sub_df = df.loc[df.index, ~df.columns.isin(used_cols)]
        sv = sub_df.mean().apply(lambda x: abs(x - self.split_value)).sort_values()
        split_id = sv.index[0]
        col_values = df.loc[df.index, split_id]
        df_index_0 = df.index[col_values == 0]
        df_index_1 = df.index[col_values == 1]
        return df_index_0, df_index_1, split_id

    def __build_tree(self, df: pd.DataFrame, used_cols: List[int], depth: int,
                     node_path: str, *args, **kwargs) -> SplitterTreeFinder.Tree:
        if len(df) < self.bucket_size or depth == 0:
            if len(df) == 0:
                Logger.log("WARNING: empty bucket")
            return self.Tree(base=self, col_id=self.leaf_id, left=None, right=None, df=df, tree_path=node_path)

        df_index_0, df_index_1, split_id = self.__split_df(df, used_cols)
        used_cols.append(split_id)
        left = self.__build_tree(df.loc[df_index_0], used_cols, depth - 1, node_path + '0', *args, **kwargs)
        right = self.__build_tree(df.loc[df_index_1], used_cols, depth - 1, node_path + '1', *args, **kwargs)
        return self.Tree(self, split_id, left, right, df, node_path, *args, **kwargs)
