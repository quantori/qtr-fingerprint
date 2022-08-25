from abc import ABC, abstractmethod
from concurrent.futures import as_completed
from typing import Generator

from substructure_finder.search_engines.search_engine import SearchEngine
from substructure_finder.search_engines.bucket_search_engine import BucketSearchEngine
from substructure_finder.fingerprint import BitFingerprint
from substructure_finder.search_engines.scope_executor import ScopeExecutor


class ExecutorSearchEngine(SearchEngine, ABC):
    @abstractmethod
    def get_executor(self) -> ScopeExecutor:
        raise NotImplementedError()

    def search(self, fingerprint: BitFingerprint) -> Generator[str, None, None]:
        executor = self.get_executor()
        buckets = self.splitter_tree.get_buckets(fingerprint)
        futures = [executor.submit(BucketSearchEngine.search_in_file, fingerprint, self.pickle_paths[bucket])
                   for bucket in buckets]
        for future in as_completed(futures):
            yield from future.result()
