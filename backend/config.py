# pylint: skip-file

import os
from dataclasses import dataclass


@dataclass
class Settings:
    APP_NAME: str = "FingerprintService"
    DEBUG_MODE = bool(os.environ["ENGINE_DEBUG"])

    INDEX_DATA_DIR = os.environ["INDEX_DATA_DIR"]
    HOST: str = os.environ["ENGINE_HOST"]
    PORT = int(os.environ["ENGINE_PORT"])


settings = Settings()
