from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor


class ScopeExecutor(Executor):
    def __del__(self):
        self.shutdown(wait=False, cancel_futures=True)


class ThreadPoolScopeExecutor(ThreadPoolExecutor, ScopeExecutor):
    pass


class ProcessPoolScopeExecutor(ProcessPoolExecutor, ScopeExecutor):
    pass
