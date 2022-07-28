import pandas as pd
from pathlib import Path
from sklearn.neighbors import BallTree
from typing import Iterable

from fp_utils.consts import PathType
from fp_utils.finders.drive_finders.drive_finder import DriveFinder
import pickle


class BallTreeFinder(DriveFinder):
    r_eps = 1e-5
    save_object = pickle

    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str = "BallTree", *args, **kwargs) -> None:
        self.tree_path = Path(directory) / (finder_name + self.file_extension)
        ball_tree = BallTree(df.values, leaf_size=2, metric='russellrao')
        self._pack((df.index, ball_tree), self.tree_path)

    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        df_index, ball_tree = self._unpack(self.tree_path)
        if not hasattr(fingerprint, 'radius'):
            fingerprint.radius = 1 - sum(fingerprint) / len(fingerprint) + self.r_eps
        ind = ball_tree.query_radius([fingerprint], fingerprint.radius)
        yield from df_index[ind[0]]
