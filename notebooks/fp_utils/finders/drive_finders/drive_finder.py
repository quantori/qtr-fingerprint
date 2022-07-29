from __future__ import annotations

from abc import ABC, abstractmethod
import pandas as pd
from pathlib import Path
from typing import Optional

from fp_utils.finders.finder import Finder
from fp_utils.consts import PathType
from fp_utils.packers import Packer


class DriveFinder(Finder, ABC):
    """Keeps data on hard drive"""

    def __new__(cls, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args,
                **kwargs) -> DriveFinder:
        """@param unique_id is used in names of files and directories associated with finder.
        NB: to simplify usage it could be skipped, though it can cause name conflicts.
        """
        super().__new__(cls)
        obj = object.__new__(cls)
        obj._unique_id = unique_id or cls.__name__
        obj._directory = Path(directory)
        obj._directory.mkdir(parents=True, exist_ok=True)
        return obj

    @abstractmethod
    def __init__(self, df: pd.DataFrame, directory: PathType, unique_id: Optional[str] = None, *args, **kwargs) -> None:
        raise NotImplementedError

    @property
    @abstractmethod
    def packer(self) -> Packer:
        raise NotImplementedError

    @property
    def data_directory(self) -> Path:
        return self._directory / self._unique_id

    @property
    def index_file(self) -> Path:
        return self._directory / (self._unique_id + self.packer.file_extension)

    def make_data_directory(self) -> None:
        self.data_directory.mkdir(parents=True, exist_ok=True)
