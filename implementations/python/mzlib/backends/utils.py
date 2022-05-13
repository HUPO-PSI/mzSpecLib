import os
import io
import gzip

from collections import deque
from typing import Any, Iterable, Optional, Union

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




class LineBuffer(object):
    lines: deque
    def __init__(self, stream: io.IOBase, lines: Iterable=None):
        if lines is None:
            lines = []
        self.lines = deque(lines)
        self.stream = stream

    def readline(self) -> bytes:
        if self.lines:
            return self.lines.popleft()
        else:
            return self.stream.readline()

    def push_line(self, line):
        self.lines.appendleft(line)

    def __iter__(self):
        while self.lines:
            yield self.lines.popleft()
        for line in self.stream:
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


def test_gzipped(f):
    """Checks the first two bytes of the
    passed file for gzip magic numbers

    Parameters
    ----------
    f : file-like or path-like

    Returns
    -------
    bool
    """
    if isinstance(f, os.PathLike):
        f = io.open(f, 'rb')
    current = f.tell()
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
    if "b" not in mode and "b" in f.mode:
        handle = io.TextIOWrapper(handle, encoding=encoding, newline=newline)
    return handle
