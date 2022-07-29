from __future__ import annotations
from abc import ABC, abstractmethod
import pandas as pd
from typing import Iterable, Optional, Generator
from fp_utils.catch_time import CatchTime


class Finder(ABC):
    def __new__(cls, *args, **kwargs) -> Finder:
        if not hasattr(cls, '__init_time_catching'):
            setattr(cls, '__init_time_catching', True)
            time_catcher = CatchTime(f'{cls.__name__} init time')
            cls.__init__ = time_catcher(cls.__init__)
        return object.__new__(cls)

    @abstractmethod
    def __init__(self, df: pd.DataFrame, *args, **kwargs) -> None:
        raise NotImplementedError

    def find(self, fingerprint: pd.Series, ans_count: Optional[int] = None) -> Generator[str, None, None]:
        if ans_count is None:
            yield from self.find_all(fingerprint)
        else:
            answer_generator = iter(self.find_all(fingerprint))
            for i in range(ans_count):
                nxt = next(answer_generator, None)
                if nxt is None:
                    break
                yield nxt

    @abstractmethod
    def find_all(self, fingerprint: pd.Series) -> Iterable[str]:
        raise NotImplementedError

    def __hash__(self):
        """Assume every finder is unique even if parameters are the same => hash == id"""
        return id(self)
