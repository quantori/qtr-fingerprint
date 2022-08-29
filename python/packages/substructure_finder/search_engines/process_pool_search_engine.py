from substructure_finder.search_engines.scope_executor import ScopeExecutor, ProcessPoolScopeExecutor
from substructure_finder.search_engines.executor_search_engine import ExecutorSearchEngine


class ProcessPoolSearchEngine(ExecutorSearchEngine):
    def get_executor(self) -> ScopeExecutor:
        return ProcessPoolScopeExecutor()
