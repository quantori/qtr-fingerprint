import sys
from concurrent.futures import as_completed
from pathlib import Path
from typing import Generator, Dict, List

from substructure_finder.search_engines.process_pool_search_engine import ProcessPoolSearchEngine
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine
from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.fingerprint import BitFingerprint


class RamProcessPoolSearchEngine(ProcessPoolSearchEngine):
    def __init__(self, db_filesystem: DbFilesystem, db_name: str):
        super().__init__(db_filesystem, db_name)
        self.bucket_search_engines: Dict[int, BucketSearchEngine] = dict()
        for bucket, file_path in self.pickle_paths.items():
            self.bucket_search_engines[bucket] = BucketSearchEngine.load(file_path)

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        executor = self.get_executor()
        buckets = list(self.splitter_tree.get_buckets(fingerprint))
        print(f'Buckets number: {len(buckets)}', file=sys.stderr)
        futures = [executor.submit(self.bucket_search_engines[bucket].search, fingerprint)
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
