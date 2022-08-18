import time
from typing import Callable

from starlette.requests import Request


async def add_process_time_header(request: Request, call_next: Callable):
    """Adds `X-Process-Time` header to the response"""
    start_time = time.perf_counter_ns()
    response = await call_next(request)
    process_time = time.perf_counter_ns() - start_time
    response.headers["X-Process-Time"] = f"{process_time / 1000000} ms."
    return response
