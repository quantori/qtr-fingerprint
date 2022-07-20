from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Any, BinaryIO, Iterable, Optional, Tuple, Generator
from sklearn.neighbors import BallTree
import joblib

from fp_utils.consts import PathType
from fp_utils.finders.drive_finder import DriveFinder
from fp_utils.logger import Logger
from fp_utils.catch_time import CatchTime


class SplitterTreeFinder(DriveFinder):
    SPLIT_VALUE = 0.5
    MAX_DEPTH = 100
    BUCKET_SIZE = 15000
    LEAF_ID = -1
    R_EPS = 1e-5
    FILE_EXTENSION = '.pickle'
    FILE_PREFIX = 'b'

    def _find(self, fingerprint: pd.Series) -> Iterable[str]:
        r = 1 - sum(fingerprint) / len(fingerprint) + self.R_EPS
        return self.__tree.find(fingerprint, r, '')
        pass

    def _load(self, file: BinaryIO) -> Any:
        return joblib.load(file)

    def _dump(self, obj: Any, file: BinaryIO) -> None:
        joblib.dump(obj, file, protocol=-1)

    def __init__(self, df: pd.DataFrame, directory: PathType) -> None:
        self.directory = Path(directory)
        self.__tree = self.__build_tree(df.index, df, list(), self.MAX_DEPTH, '')

    def _get_file_path(self, node_path: str) -> Path:
        return self.directory / (self.FILE_PREFIX + node_path + self.FILE_EXTENSION)

    # @CatchTime('pack bucket')
    def _pack_bucket(self, df_index: pd.Index, df: pd.DataFrame, node_path: str) -> None:
        file_path = self._get_file_path(node_path)
        index = df_index, BallTree(df.loc[df_index].values, leaf_size=2, metric='russellrao')
        self._pack(index, file_path)

    # @CatchTime('unpack bucket')
    def _unpack_bucket(self, node_path: str) -> Tuple[pd.Index, BallTree]:
        file_path = self._get_file_path(node_path)
        return self._unpack(file_path)

    class Tree:
        def __init__(self, base: SplitterTreeFinder, col_id: int, left: Optional[SplitterTreeFinder.Tree],
                     right: Optional[SplitterTreeFinder.Tree], df_index: pd.Index, df: pd.DataFrame,
                     node_path: str) -> None:
            self.base = base
            self.col_id = col_id
            self.left = left
            self.right = right
            if self.is_leaf:
                self.base._pack_bucket(df_index, df, node_path)

        @property
        def is_leaf(self) -> bool:
            return self.col_id == self.base.LEAF_ID

        def find(self, fingerprint: pd.Series, r: int, node_path: str) -> Generator[str, None, None]:
            if self.is_leaf:
                index, ball_tree = self.base._unpack_bucket(node_path)
                ind = ball_tree.query_radius([fingerprint], r)
                yield from index[ind[0]]
            else:
                if fingerprint[self.col_id] == 0:
                    yield from self.left.find(fingerprint, r, node_path + '0')
                yield from self.right.find(fingerprint, r, node_path + '1')

    # @CatchTime("split df")
    def __split_df(self, df_index: pd.Index, df: pd.DataFrame, used_cols: List[int]) -> Tuple[pd.Index, pd.Index, int]:
        sub_df = df.loc[df_index, ~df.columns.isin(used_cols)]
        sv = sub_df.mean().apply(lambda x: abs(x - self.SPLIT_VALUE)).sort_values()
        split_id = sv.index[0]
        col_values = df.loc[df_index, split_id]
        df_index_0 = df_index[col_values == 0]
        df_index_1 = df_index[col_values == 1]
        return df_index_0, df_index_1, split_id

    def __build_tree(self, df_index: pd.Index, df: pd.DataFrame, used_cols: List[int], depth: int,
                     node_path: str) -> SplitterTreeFinder.Tree:
        if df_index.size < self.BUCKET_SIZE or depth == 0:
            if df_index.size == 0:
                Logger.log("WARNING: empty bucket")
            return self.Tree(base=self, col_id=self.LEAF_ID, left=None, right=None, df_index=df_index, df=df,
                             node_path=node_path)

        df_index_0, df_index_1, split_id = self.__split_df(df_index, df, used_cols)
        used_cols.append(split_id)
        left = self.__build_tree(df_index_0, df, used_cols, depth - 1, node_path + '0')
        right = self.__build_tree(df_index_1, df, used_cols, depth - 1, node_path + '1')
        return self.Tree(base=self, col_id=split_id, left=left, right=right, df_index=df_index, df=df,
                         node_path=node_path)
