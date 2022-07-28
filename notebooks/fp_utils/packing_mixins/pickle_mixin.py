import pickle
from pathlib import Path
from typing import Any

from fp_utils.consts import PathType
from fp_utils.packing_mixins.packing_mixin import PackingMixin


class PickleMixin(PackingMixin):
    file_extension = '.pickle'

    @staticmethod
    def _pack(obj: object, file_path: PathType) -> None:
        with Path(file_path).open('wb') as f:
            pickle.dump(obj, f, protocol=-1)

    @staticmethod
    def _unpack(file_path: PathType) -> Any:
        with Path(file_path).open('rb') as f:
            obj = pickle.load(f)
        return obj
