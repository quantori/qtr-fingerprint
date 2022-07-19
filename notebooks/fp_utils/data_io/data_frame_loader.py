import numpy as np
import pandas as pd

from typing import Union
from fp_utils.consts import PathType


class DataFrameLoader:
    @staticmethod
    def pickle(df_path: str, columns: Union[np.ndarray, PathType, None] = None) -> pd.DataFrame:
        df = pd.read_pickle(df_path)
        return df if columns is None else DataFrameLoader.__pick_columns(df, columns)

    @staticmethod
    def __pick_columns(df: pd.DataFrame, columns: Union[np.ndarray, PathType]) -> pd.DataFrame:
        if isinstance(columns, PathType):
            columns = np.load(str(columns))
        assert isinstance(columns, np.ndarray)
        return df[columns]
