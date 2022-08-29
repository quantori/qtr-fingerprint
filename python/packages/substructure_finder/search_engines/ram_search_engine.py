import sys
from concurrent.futures import as_completed, ThreadPoolExecutor
from typing import Generator, Dict, List

from substructure_finder.search_engines.search_engine import SearchEngine
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine
from substructure_finder.db_filesystem import DbFilesystem
from substructure_finder.fingerprint import BitFingerprint


class RamSearchEngine(SearchEngine):
    def __init__(self, db_filesystem: DbFilesystem, db_name: str):
        super().__init__(db_filesystem, db_name)
        self.bucket_search_engines: Dict[int, BucketSearchEngine] = dict()
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(self._load_buckets, buckets) for buckets in
                       db_filesystem.base_dir_to_buckets(db_name).values()]
            for future in as_completed(futures):
                self.bucket_search_engines.update(future.result())

    def _load_buckets(self, buckets: List[int]) -> Dict[int, BucketSearchEngine]:
        return dict((bucket, BucketSearchEngine.load(self.pickle_paths[bucket])) for bucket in buckets)

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        buckets = list(self.splitter_tree.get_buckets(fingerprint))
        print(f'Buckets number: {len(buckets)}', file=sys.stderr)
        for bucket in buckets:
            yield from self.bucket_search_engines[bucket].search(fingerprint)
