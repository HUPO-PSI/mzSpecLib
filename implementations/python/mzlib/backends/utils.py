import os
import io
import gzip

from collections import deque
from urllib import parse as urlparse
from typing import Any, Dict, Iterable, Mapping, Optional, Union

DEFAULT_BUFFER_SIZE = int(2e6)
GZIP_MAGIC = b'\037\213'

GzipFile = gzip.GzipFile

try:
    # Fast random acces with Gzip compatibility
    import idzip
    idzip.compressor.IdzipWriter.enforce_extension = False

    GzipFile = idzip.IdzipFile
except ImportError:
    pass


class _LineBuffer(object):
    """
    An implementation detail that treats a stream/iterator over line strings as LIFO
    queue that can have lines pushed back onto it.
    """

    lines: deque
    stream: io.IOBase
    last_line: str
    _stream_is_file_like: bool

    def __init__(self, stream: io.IOBase, lines: Iterable=None, last_line: str=None):
        if lines is None:
            lines = []
        self.lines = deque(lines)
        self.stream = stream
        self.last_line = last_line
        self._stream_is_file_like = hasattr(self.stream, 'readline')

    def readline(self) -> Union[bytes, str]:
        if self.lines:
            line = self.lines.popleft()
        else:
            line = self.stream.readline() if self._stream_is_file_like else next(self.stream)
        self.last_line = line
        return line

    def push_line(self, line=None):
        if line is None:
            line = self.last_line
            self.last_line = None
        if line is None:
            raise ValueError("Cannot push empty value after the backtrack line is consumed")
        self.lines.appendleft(line)

    def __iter__(self):
        while self.lines:
            line = self.lines.popleft()
            self.last_line = line
            yield line
        for line in self.stream:
            self.last_line = line
            yield line

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def try_cast(value: Any) -> Union[str, int, float, Any]:
    if value is None:
        return value
    if not isinstance(value, str):
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def test_gzipped(f) -> bool:
    """
    Checks the first two bytes of the
    passed file for gzip magic numbers

    Parameters
    ----------
    f : file-like or path-like
        The file to test

    Returns
    -------
    bool
    """
    if isinstance(f, os.PathLike):
        f = io.open(f, 'rb')
    try:
        current = f.tell()
    except OSError:
        return False
    f.seek(0)
    magic = f.read(2)
    f.seek(current)
    return magic == GZIP_MAGIC


def starts_with_gz_magic(bytestring):
    '''Tests whether or not a byte string starts with
    the GZIP magic bytes.

    Parameters
    ----------
    bytestring : bytes
        The bytes to test.

    Returns
    -------
    bool
    '''
    return bytestring.startswith(GZIP_MAGIC)


def open_stream(f: Union[io.IOBase, os.PathLike], mode='rt', buffer_size: Optional[int]=None, encoding: Optional[str]='utf8', newline=None):
    '''Select the file reading type for the given path or stream.

    Detects whether the file is gzip encoded.
    '''
    if buffer_size is None:
        buffer_size = DEFAULT_BUFFER_SIZE
    if 'r' in mode:
        if not hasattr(f, 'read'):
            f = io.open(f, 'rb')
        # On Py2, dill doesn't behave correctly with io-derived objects, so we have to
        # patch it below. Don't try to wrap an io.TextIOWrapper on Py3.
        if not isinstance(f, io.BufferedReader) and not isinstance(f, io.TextIOWrapper):
            buffered_reader = io.BufferedReader(f, buffer_size)
        else:
            buffered_reader = f
        if test_gzipped(buffered_reader):
            handle = GzipFile(fileobj=buffered_reader, mode='rb')
        else:
            handle = buffered_reader
    else:
        raise NotImplementedError("Haven't implemented automatic output stream determination")
    try:
        fmode = f.mode
    except AttributeError:
        fmode = 'b'
    if "b" not in mode and "b" in fmode:
        handle = io.TextIOWrapper(handle, encoding=encoding, newline=newline)
    return handle


class CaseInsensitiveDict(Dict[str, Any]):
    """A case sensitive version of a dictionary with string keys."""

    def __init__(self, base=None, **kwargs):
        if base is not None:
            self.update(base)
        if kwargs:
            self.update(kwargs)

    def __getitem__(self, key: str):
        return super().__getitem__(key.lower())

    def __setitem__(self, key: str, value):
        super().__setitem__(key.lower(), value)

    def __delitem__(self, key: str):
        super().__delitem__(key.lower())

    def __contains__(self, __o: str) -> bool:
        return super().__contains__(__o.lower())

    def get(self, key: str, default=None):
        return super().get(key.lower(), default)

    def update(self, value: Mapping[str, Any]):
        super().update({k.lower(): v for k, v in value.items()})


def urlify(path: str) -> str:
    """Convert a path into a URL if it is not already one."""
    parsed = urlparse.urlparse(path)
    if parsed.scheme == '':
        parsed = parsed._replace(scheme='file')
    return urlparse.urlunparse(parsed)