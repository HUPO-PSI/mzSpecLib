import warnings
import logging

from typing import Any, Dict, Iterator, Optional, List, DefaultDict, Union

from numbers import Integral
from collections import defaultdict

from mzlib.index.base import IndexRecordBase

from .base import IndexBase, IndexRecordBase

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class _IndexAttr:
    __slots__ = ()

    def get(self, key: str, default=None) -> Any:
        if self.attributes is not None:
            return self.attributes.get(key, default)
        return default

    def set(self, key: str, value: Any):
        if self.attributes is not None:
            self.attributes[key] = value
        else:
            self.attributes = {key: value}


class IndexRecord(IndexRecordBase, _IndexAttr):
    """
    A spectrum index record.

    Attributes
    ----------
    number : int
        A numerical identifier for the spectrum
    offset : int
        The offset in the file to reach the spectrum (in bytes if appropriate)
    name : str,
        A text identifier for this spectrum.
    analyte : str, optional
        A text representation of the analyte for that record
    attributes : Dict[str, Any], optional
        A key-value pair collection of this record.
    """

    __slots__ = ('number', 'offset', 'name', 'analyte', 'index', 'attributes')

    number: int
    offset: int
    name: str
    index: int
    analyte: Any
    attributes: Optional[Dict[str, Any]]

    def __init__(self, number, offset, name, analyte, index: int=None, attributes=None):
        self.number = number
        self.offset = offset
        self.name = name
        self.analyte = analyte
        self.index = index
        self.attributes = attributes

    def __repr__(self):
        template = f"{self.__class__.__name__}({self.number}, {self.offset}, {self.name}, {self.analyte}, {self.index}, {self.attributes})"
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

    def to_dict(self) -> Dict:
        return {
            k: getattr(self, k, None) for k in self.__slots__
        }

    @classmethod
    def from_dict(cls, state: Dict) -> 'IndexRecord':
        return cls(**state)


class ClusterIndexRecord(IndexRecordBase, _IndexAttr):
    """
    A spectrum cluster index record.

    Attributes
    ----------
    number : int
        A numerical identifier for this spectrum.
    offset : int
        The offset in the file to reach the spectrum (in bytes if appropriate)
    attributes : Dict[str, Any], optional
        A key-value pair collection of this record
    """

    __slots__ = ('number', 'offset', 'attributes')

    def __init__(self, number, offset, attributes=None):
        self.number = number
        self.offset = offset
        self.attributes = attributes

    def __repr__(self):
        template = f"{self.__class__.__name__}({self.number}, {self.offset}, {self.attributes})"
        return template

    def __eq__(self, other):
        if self.number != other.number:
            return False
        elif self.offset != other.offset:
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

    def to_dict(self) -> Dict:
        return {
            k: getattr(self, k, None) for k in self.__slots__
        }

    @classmethod
    def from_dict(cls, state: Dict) -> 'ClusterIndexRecord':
        return cls(**state)


class MemoryIndex(IndexBase):
    records: List[IndexRecord]
    cluster_records: List[ClusterIndexRecord]
    metadata: Dict[str, Any]

    _dirty: bool
    _by_key: Dict[int, IndexRecord]
    _by_name: DefaultDict[str, List[IndexRecord]]
    _by_attr: DefaultDict[str, DefaultDict[Any, List[IndexRecord]]]

    @classmethod
    def from_filename(cls, filename, library=None):
        inst = cls()
        return inst, False

    def __init__(self, records=None, cluster_records=None, metadata=None):
        self.records = list(records or [])
        self.cluster_records = list(cluster_records or [])
        self._by_name = defaultdict(list)
        self._by_key = {}
        self._by_attr = defaultdict(lambda: defaultdict(list))
        self.metadata = metadata or {}
        self._dirty = True

    def iter_clusters(self) -> Iterator[IndexRecordBase]:
        """Iterate over cluster entries in the index."""
        return iter(self.cluster_records)

    def iter_spectra(self):
        """Iterate over spectrum entries in the index."""
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
                return self._by_key[i]
            except IndexError as err:
                raise KeyError(i) from err
        elif isinstance(i, slice):
            start = i.start
            stop = i.stop
            if start is None:
                start = min(self._by_key) if self._by_key else 0
            if stop is None:
                stop = max(self._by_key) if self._by_key else 0
            return [self._by_key[i] for i in range(start, stop) if i in self._by_key]
        if i in self._by_name:
            records = self._by_name[i]
            if len(records) == 1:
                return records[0]
            else:
                return records
        else:
            raise KeyError(i)

    def search_clusters(self, i=None, **kwargs):
        if self._dirty:
            self._update_index()
        if i is None and kwargs:
            # Executing attribute query
            raise NotImplementedError()
        if isinstance(i, Integral):
            try:
                return self.cluster_records[i]
            except IndexError as err:
                raise KeyError(i) from err
        elif isinstance(i, slice):
            return self.cluster_records[i]

    def __getitem__(self, i):
        return self._get_by_index(i)

    def _get_by_index(self, i: Union[int, slice]) -> Union[IndexRecord, List[IndexRecord]]:
        return self.records[i]

    def _update_index(self):
        self.records.sort(key=lambda x: x.number)

        self._by_name = defaultdict(list)
        for record in self:
            self._by_key[record.number] = record
            self._by_name[record.name].append(record)
        self._dirty = False

    def add(self, number: int, offset: int, name: str, analyte: Any, attributes=None):
        """
        Add a new entry to the spectrum index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        name : str,
            A text identifier for this spectrum.
        analyte : str, optional
            A text representation of the analyte for that record
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        n = len(self.records)
        record = IndexRecord(number, offset, name, analyte, n, attributes)
        self.records.append(record)
        self._dirty = True

    def add_cluster(self, number: int, offset: int, attributes=None):
        """
        Add a new entry to the spectrum index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        record = ClusterIndexRecord(number, offset, attributes)
        self.cluster_records.append(record)
        self._dirty = True


    def commit(self):
        self._update_index()
