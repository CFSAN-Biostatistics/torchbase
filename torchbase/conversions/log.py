"""Shared logging setup for torchbase converters."""

import logging
import sys

TRACE = 5
logging.addLevelName(TRACE, "TRACE")


class _ConverterFormatter(logging.Formatter):
    """Indent output visually by severity to show hierarchy."""

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno >= logging.INFO:
            prefix = ""
        elif record.levelno >= logging.DEBUG:
            prefix = "  · "
        else:  # TRACE
            prefix = "    · "
        record.msg = prefix + str(record.msg)
        return super().format(record)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"torchbase.conversions.{name}")


def setup_logging(verbosity: int) -> None:
    """Configure the torchbase logger for the given verbosity count.

    0 → silent (WARNING, nothing from conversions)
    1 → INFO  (-v)
    2 → DEBUG (-vv)
    3+ → TRACE (-vvv)
    """
    if verbosity == 0:
        return
    level = logging.INFO if verbosity == 1 else (logging.DEBUG if verbosity == 2 else TRACE)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(_ConverterFormatter("%(message)s"))

    logger = logging.getLogger("torchbase")
    logger.setLevel(level)
    if not any(isinstance(h, logging.StreamHandler) and h.stream is sys.stdout
               for h in logger.handlers):
        logger.addHandler(handler)
    logger.propagate = False
