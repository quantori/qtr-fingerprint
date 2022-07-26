from __future__ import annotations

from abc import ABC, abstractmethod
import pandas as pd
from typing import Union, Any, BinaryIO
from pathlib import Path
import pickle

from fp_utils.finders.finder import Finder
from fp_utils.consts import PathType


class DriveFinder(Finder, ABC):
    """Keeps data on hard drive"""

    save_object = pickle

    @property
    def file_extension(self):
        return f'.{self.save_object.__name__}'

    def __new__(cls, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> DriveFinder:
        super().__new__(cls)
        Path(directory).mkdir(parents=True, exist_ok=True)
        return object.__new__(cls)

    @abstractmethod
    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        raise NotImplementedError

    def _pack(self, obj: Any, file: PathType) -> None:
        with Path(file).open('wb') as f:
            self._dump(obj, f)

    def _unpack(self, file: Union[str, Path]) -> Any:
        with Path(file).open('rb') as f:
            obj = self._load(f)
        return obj

    def _dump(self, obj: Any, file: BinaryIO) -> None:
        self.save_object.dump(obj, file, protocol=-1)

    def _load(self, file: BinaryIO) -> Any:
        return self.save_object.load(file)
