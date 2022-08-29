from substructure_finder.search_engines.scope_executor import ScopeExecutor, ThreadPoolScopeExecutor
from substructure_finder.search_engines.executor_search_engine import ExecutorSearchEngine


class ThreadPoolSearchEngine(ExecutorSearchEngine):
    def get_executor(self) -> ScopeExecutor:
        return ThreadPoolScopeExecutor()
