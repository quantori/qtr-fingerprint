from pathlib import Path
from typing import Literal, Iterable, Dict
from enum import IntEnum
from concurrent.futures import as_completed, ProcessPoolExecutor

from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.molecule import Molecule
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine


class DbType(IntEnum):
    RawDB = 0
    DB = 1


class DbStat:
    def __init__(self, db_filesystem: DbFilesystem, db_name: str, db_type: Literal["raw_db", "db"]):
        self.db_filesystem = db_filesystem
        self.db_name = db_name
        self.db_type = DbType.DB if db_type == 'db' else DbType.RawDB

    @property
    def size(self) -> int:
        if self.db_type == DbType.DB:
            file_paths = self.db_filesystem.pickle_paths(self.db_name)
            count_func = DbStat._count_molecules_in_pickles
        elif self.db_type == DbType.RawDB:
            file_paths = self.db_filesystem.rb_paths(self.db_name)
            count_func = DbStat._count_molecules_in_rbs
        else:
            assert False, "Wrong db_type"

        molecules_number = 0
        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(count_func, file_paths, buckets) for buckets in
                       self.db_filesystem.base_dir_to_buckets(self.db_name).values()]
            for future in as_completed(futures):
                molecules_number += future.result()
        return molecules_number

    @staticmethod
    def _count_molecules_in_pickles(file_paths: Dict[int, Path], buckets: Iterable[int]) -> int:
        molecules_number = 0
        for bucket in buckets:
            bucket_search_engine = BucketSearchEngine.load(file_paths[bucket])
            molecules_number += len(bucket_search_engine.molecules)
        return molecules_number

    @staticmethod
    def _count_molecules_in_rbs(file_paths: Dict[int, Path], buckets: Iterable[int]) -> int:
        molecules_number = 0
        for bucket in buckets:
            with file_paths[bucket].open('rb') as stream:
                molecules = Molecule.load_molecules_from_rb_stream(stream)
            molecules_number += len(molecules)
        return molecules_number
