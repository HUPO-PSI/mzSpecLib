import logging
import re
from typing import Dict


class LevelAwareColoredLogFormatter(logging.Formatter):
    try:
        from colorama import Fore, Style
        # GREY = Fore.WHITE
        GREY = ''
        BLUE = Fore.BLUE
        GREEN = Fore.GREEN
        YELLOW = Fore.YELLOW
        RED = Fore.RED
        BRIGHT = Style.BRIGHT
        DIM = Style.DIM
        BOLD_RED = Fore.RED + Style.BRIGHT
        RESET = Style.RESET_ALL
    except ImportError:
        GREY = ''
        BLUE = ''
        GREEN = ''
        YELLOW = ''
        RED = ''
        BRIGHT = ''
        DIM = ''
        BOLD_RED = ''
        RESET = ''

    def _colorize_field(self, fmt: str, field: str, color: str) -> str:
        return re.sub("(" + field + ")", color + r"\1" + self.RESET, fmt)

    def _patch_fmt(self, fmt: str, level_color: str) -> str:
        fmt = self._colorize_field(fmt, r"%\(asctime\)s", self.GREEN)
        fmt = self._colorize_field(fmt, r"%\(name\).*?s", self.BLUE)
        fmt = self._colorize_field(fmt, r"%\(message\).*?s", self.GREY)
        if level_color:
            fmt = self._colorize_field(fmt, r"%\(levelname\).*?s", level_color)
        return fmt

    def __init__(self, fmt, level_color=None, **kwargs):
        fmt = self._patch_fmt(fmt, level_color=level_color)
        super().__init__(fmt, **kwargs)


class ColoringFormatter(logging.Formatter):
    level_to_color = {
        logging.INFO: LevelAwareColoredLogFormatter.GREEN,
        logging.DEBUG: LevelAwareColoredLogFormatter.GREY + LevelAwareColoredLogFormatter.DIM,
        logging.WARN: LevelAwareColoredLogFormatter.YELLOW + LevelAwareColoredLogFormatter.BRIGHT,
        logging.ERROR: LevelAwareColoredLogFormatter.BOLD_RED,
        logging.CRITICAL: LevelAwareColoredLogFormatter.BOLD_RED,
        logging.FATAL: LevelAwareColoredLogFormatter.RED + LevelAwareColoredLogFormatter.DIM,
    }

    _formatters: Dict[int, LevelAwareColoredLogFormatter]

    def __init__(self, fmt: str, **kwargs):
        self._formatters = {}
        for level, style in self.level_to_color.items():
            self._formatters[level] = LevelAwareColoredLogFormatter(fmt, level_color=style, **kwargs)

    def format(self, record: logging.LogRecord) -> str:
        fmtr = self._formatters[record.levelno]
        return fmtr.format(record)
