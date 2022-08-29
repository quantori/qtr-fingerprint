from asyncio import as_completed
from typing import Generator

from substructure_finder.search_engines.process_pool_search_engine import ProcessPoolSearchEngine
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine
from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.fingerprint import BitFingerprint


class RamProcessPoolSearchEngine(ProcessPoolSearchEngine):
    def __init__(self, db_filesystem: DbFilesystem, db_name: str):
        super().__init__(db_filesystem, db_name)
        self.bucket_search_engines = dict((bucket, BucketSearchEngine.load(file_path))
                                          for bucket, file_path in self.pickle_paths.items())

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        executor = self.get_executor()
        buckets = self.splitter_tree.get_buckets(fingerprint)
        futures = [executor.submit(BucketSearchEngine.search, self.bucket_search_engines[bucket], fingerprint)
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
