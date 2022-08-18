import logging

from api.routers import router as api_router
from config import settings
from config.middlewares import add_process_time_header

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from starlette.middleware.base import BaseHTTPMiddleware
from uvicorn import run

log = logging.getLogger(__name__)

app = FastAPI(
    docs_url="/api/docs",
    openapi_url="/api/openapi.json",
)

app.add_middleware(
    CORSMiddleware,
    allow_methods=["POST", "GET", "OPTIONS", "DELETE", "PUT"],
    allow_headers=["*"],
)

if settings.DEBUG_MODE:
    # Adding some parameters for debugging purpose
    app.add_middleware(BaseHTTPMiddleware, dispatch=add_process_time_header)

app.include_router(api_router, tags=["api"], prefix="/api")

if __name__ == "__main__":
    run(
        "main:app",
        host=settings.HOST,
        reload=settings.DEBUG_MODE,
        port=settings.PORT,
        log_config=settings.LOG_CONFIG,
    )
