from abc import ABC, abstractmethod
from typing import Any

from fp_utils.consts import PathType


class Packer(ABC):
    @staticmethod
    @abstractmethod
    def pack(obj: object, file_path: PathType) -> None:
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def unpack(file_path: PathType) -> Any:
        raise NotImplementedError

    @property
    @abstractmethod
    def file_extension(self):
        raise NotImplementedError
