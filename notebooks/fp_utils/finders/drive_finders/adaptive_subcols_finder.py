from abc import abstractmethod, ABC
import pandas as pd
from typing import Iterable

from fp_utils.consts import PathType
from fp_utils.finders.drive_finders.subcols_finder import SubColsFinder


class AdaptiveSubColsFinder(SubColsFinder, ABC):
    @staticmethod
    @abstractmethod
    def columns_chooser(df: pd.DataFrame) -> Iterable[int]:
        raise NotImplementedError

    @property
    def columns(self) -> Iterable[int]:
        return self._columns

    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        self._columns = self.__class__.columns_chooser(df)
        super().__init__(df, directory, finder_name, *args,  **kwargs)
