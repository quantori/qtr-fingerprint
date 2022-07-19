from fp_utils.finders.finder import Finder
from abc import ABC, abstractmethod
import pandas as pd
from typing import Union, Any, BinaryIO
from pathlib import Path
from fp_utils.consts import PathType


class DriveFinder(Finder, ABC):
    """Keeps data on hard drive"""

    @abstractmethod
    def __init__(self, df: pd.DataFrame, directory: PathType, *args, **kwargs) -> None:
        raise NotImplementedError

    def _pack(self, obj: Any, file: PathType) -> None:
        with Path(file).open('wb') as f:
            self._dump(obj, f)

    def _unpack(self, file: Union[str, Path]) -> Any:
        with Path(file).open('rb') as f:
            obj = self._load(f)
        return obj

    @abstractmethod
    def _dump(self, obj: Any, file: BinaryIO) -> None:
        raise NotImplementedError

    @abstractmethod
    def _load(self, file: BinaryIO) -> Any:
        raise NotImplementedError
