import warnings
import logging

from numbers import Integral
from collections import defaultdict

from .base import IndexBase

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IndexRecord(object):
    __slots__ = ('number', 'offset', 'name', 'analyte', 'attributes')

    def __init__(self, number, offset, name, analyte, attributes=None):
        self.number = number
        self.offset = offset
        self.name = name
        self.analyte = analyte
        self.attributes = attributes or {}

    def __repr__(self):
        template = f"{self.__class__.__name__}({self.number}, {self.offset}, {self.name}, {self.analyte}, {self.attributes})"
        return template

    def __eq__(self, other):
        if self.number != other.number:
            return False
        elif self.offset != other.offset:
            return False
        elif self.name != other.name:
            return False
        elif self.analyte != other.analyte:
            return False
        if bool(self.attributes) == bool(other.attributes):
            if bool(self.attributes) and self.attributes != other.attributes:
                return False
            # Implicitly allow None and empty dictionaries to be the same
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.name)


class MemoryIndex(IndexBase):

    @classmethod
    def from_filename(cls, filename, library=None):
        inst = cls()
        return inst, False

    def __init__(self, records=None, metadata=None):
        self.records = list(records or [])
        self._by_name = defaultdict(list)
        self._by_attr = defaultdict(lambda: defaultdict(list))
        self.metadata = metadata or {}
        self._dirty = True

    def __iter__(self):
        return iter(self.records)

    def __len__(self):
        return len(self.records)

    def search(self, i=None, **kwargs):
        if self._dirty:
            self._update_index()
        if i is None and kwargs:
            # Executing attribute query
            raise NotImplementedError()
        if isinstance(i, Integral):
            try:
                return self.records[i]
            except IndexError as err:
                raise KeyError(i) from err
        elif isinstance(i, slice):
            return self.records[i]
        if i in self._by_name:
            records = self._by_name[i]
            if len(records) == 1:
                return records[0]
            else:
                return records
        else:
            raise KeyError(i)

    def __getitem__(self, i):
        return self.search(i)

    def _update_index(self):
        self.records.sort(key=lambda x: x.number)

        self._by_name = defaultdict(list)
        for record in self:
            self._by_name[record.name].append(record)

        self._dirty = False

    def add(self, number, offset, name, analyte, attributes=None):
        record = IndexRecord(number, offset, name, analyte, attributes)
        self.records.append(record)
        self._dirty = True

    def commit(self):
        self._update_index()
