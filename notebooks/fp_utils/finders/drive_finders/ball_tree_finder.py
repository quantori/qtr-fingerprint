from abc import ABC

import pandas as pd
from sklearn.neighbors import BallTree
from typing import Iterable, Optional

from fp_utils.consts import PathType
from fp_utils.finders.drive_finders.drive_finder import DriveFinder


class BallTreeFinder(DriveFinder, ABC):
    r_eps = 1e-5

    def __init__(self, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args, **kwargs) -> None:
        ball_tree = BallTree(df.values, leaf_size=2, metric='russellrao')
        self.packer.pack((df.index, ball_tree), self.index_file)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df_index, ball_tree = self.packer.unpack(self.index_file)
        if not hasattr(fingerprint, 'radius'):
            fingerprint.radius = 1 - sum(fingerprint) / len(fingerprint) + self.r_eps
        ind = ball_tree.query_radius([fingerprint], fingerprint.radius)
        yield from df_index[ind[0]]
