import sys


class Logger:
    @staticmethod
    def log(text: object, verbose: bool = True, file=sys.stdout) -> None:
        if not verbose:
            return
        print(text, file=file)
