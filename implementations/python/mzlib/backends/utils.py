from collections import deque

class LineBuffer(object):
    def __init__(self, stream, lines=None):
        if lines is None:
            lines = []
        self.lines = deque(lines)
        self.stream = stream

    def readline(self):
        if self.lines:
            self.lines.popleft()
        else:
            return self.stream.readline()

    def push_line(self, line):
        self.lines.appendleft(line)

    def __iter__(self):
        while self.lines:
            yield self.lines.popleft
        for line in self.stream:
            yield line

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


def try_cast(value):
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