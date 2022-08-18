# pylint: skip-file

import os

from pydantic import BaseSettings


class Settings(BaseSettings):
    APP_NAME: str = "FingerprintService"
    PROTOCOL: str = "http"
    HOST: str = os.environ["ENGINE_HOST"]
    PORT: int = int(os.environ["ENGINE_PORT"])
    URL: str = f"{PROTOCOL}://{HOST}:{PORT}"
    DEBUG_MODE: bool = bool(os.environ["DEBUG_MODE"])
    LOG_CONFIG: str = "log_config.yaml"

    INDEX_DATA_DIR = os.environ["INDEX_DATA_DIR"]
    CORS_ORIGINS_ALLOWED: list[str] = ["http://localhost:8000"]


settings = Settings()
