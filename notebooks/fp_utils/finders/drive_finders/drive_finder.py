from __future__ import annotations

from abc import ABC, abstractmethod
import pandas as pd
from pathlib import Path

from fp_utils.finders.finder import Finder
from fp_utils.consts import PathType
from fp_utils.packing_mixins.packing_mixin import PackingMixin


class DriveFinder(Finder, PackingMixin, ABC):
    """Keeps data on hard drive"""

    def __new__(cls, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> DriveFinder:
        super().__new__(cls)
        Path(directory).mkdir(parents=True, exist_ok=True)
        return object.__new__(cls)

    @abstractmethod
    def __init__(self, df: pd.DataFrame, directory: PathType, finder_name: str, *args, **kwargs) -> None:
        raise NotImplementedError
