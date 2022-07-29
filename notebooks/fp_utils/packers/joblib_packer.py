from pathlib import Path
from typing import Any
import joblib

from fp_utils.consts import PathType
from fp_utils.packers.packer import Packer


class JoblibPacker(Packer):
    file_extension = '.joblib'

    @staticmethod
    def pack(obj: object, file_path: PathType) -> None:
        with Path(file_path).open('wb') as f:
            joblib.dump(obj, f, protocol=-1)

    @staticmethod
    def unpack(file_path: PathType) -> Any:
        with Path(file_path).open('rb') as f:
            obj = joblib.load(f)
        return obj
