import os
from pathlib import Path
from typing import List, Dict


class DbFilesystem:
    def __init__(self, data_directories: List[Path], other_directory: Path):
        self._data_directories = [Path(directory) for directory in data_directories]
        self._other_directory = Path(other_directory)

    def bucket_dir_paths(self, db_name: str) -> Dict[int, Path]:
        result = dict()
        for data_dir in self._data_directories:
            db_dir = data_dir / db_name
            for name in map(str, os.listdir(db_dir)):
                result[int(name)] = db_dir / name
        return result

    def base_dir_to_buckets(self, db_name: str) -> Dict[Path, List[int]]:
        result = dict()
        for bucket, bucket_path in self.bucket_dir_paths(db_name).items():
            result.setdefault(self.bucket_base_dir(bucket_path.absolute()), []).append(bucket)
        return result

    def tree_path(self, db_name: str) -> Path:
        return self._other_directory / db_name / 'tree'

    def pickle_paths(self, db_name: str) -> Dict[int, Path]:
        bucket_dirs = self.bucket_dir_paths(db_name)
        return dict((bucket, self.pickle_file_from_bucket_dir(directory)) for bucket, directory in bucket_dirs.items())

    @property
    def data_directories(self):
        return self._data_directories

    @property
    def other_directory(self):
        return self._other_directory

    @classmethod
    def bucket_base_dir(cls, bucket_path: Path) -> Path:
        return bucket_path.absolute().parent.parent

    @classmethod
    def rb_file_from_bucket_dir(cls, bucket_dir: Path) -> Path:
        return bucket_dir / 'data.rb'

    @classmethod
    def col_file_from_bucket_dir(cls, bucket_dir: Path) -> Path:
        return bucket_dir / 'columns.col'

    @classmethod
    def pickle_file_from_bucket_dir(cls, bucket_dir: Path) -> Path:
        return bucket_dir / 'data.pickle'

    def delete_db(self, db_name: str):
        del_dirs = [directory / db_name for directory in self._data_directories]
        del_dirs += [self.other_directory / db_name]
        print('Do you want to delete this directories:\n' + '\n'.join(map(str, del_dirs)))
        user_answer = input('Y/n?: ')
        if user_answer.lower() != 'y':
            print('abort operation')
            return
        for directory in del_dirs:
            os.system(f"rm -r {directory}")
