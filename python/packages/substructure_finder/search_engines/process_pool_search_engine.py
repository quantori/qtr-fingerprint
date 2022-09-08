from substructure_finder.scope_executor import ScopeExecutor, ProcessPoolScopeExecutor
from substructure_finder.search_engines.executor_search_engine import ExecutorSearchEngine


class ProcessPoolSearchEngine(ExecutorSearchEngine):
    def get_executor(self) -> ProcessPoolScopeExecutor:
        return ProcessPoolScopeExecutor()
