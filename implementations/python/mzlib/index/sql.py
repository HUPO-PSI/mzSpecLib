import io
import os
import numbers
import pathlib
import logging
from typing import Iterator, List, Union

from sqlalchemy import Column, ForeignKey, Integer, Float, String, DateTime, Text, LargeBinary

from mzlib.index.base import IndexRecordBase
try: # For SQLAlchemy 2.0
    from sqlalchemy.orm import declarative_base
except ImportError:
    from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.engine import Engine
from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker, scoped_session

from .base import IndexBase


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


Base = declarative_base()

#### Define the database tables as classes


class SpectrumLibraryIndexAttribute(Base):
    __tablename__ = 'spectrum_library_index_attribute'
    id = Column(Integer, primary_key=True)
    name = Column(String(255), nullable=False)
    value = Column(String(1024), nullable=False)


#### Define the database tables as classes
class SpectrumLibraryIndexRecord(Base):
    __tablename__ = 'spectrum_library_index_record'
    id = Column(Integer, primary_key=True)
    number = Column(Integer, nullable=False, index=True)
    offset = Column(Integer, nullable=False)
    name = Column(String(1024), nullable=False)
    index = Column(Integer, nullable=False, index=True)
    analyte = Column(String(2014), nullable=True)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.number}, {self.offset}, {self.name}, {self.analyte})"


class ClusterSpectrumLibraryIndexRecord(Base):
    __tablename__ = 'cluster_spectrum_library_index_record'
    id = Column(Integer, primary_key=True)
    number = Column(Integer, nullable=False, index=True)
    offset = Column(Integer, nullable=False)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.number}, {self.offset}, {self.name}, {self.analyte})"


class SQLIndex(IndexBase):
    extension = '.splindex'

    filename: str
    index_filename: str
    _cache: SpectrumLibraryIndexRecord
    session: scoped_session
    engine: Engine

    @classmethod
    def from_filename(cls, filename, library=None):
        if not isinstance(filename, (str, pathlib.Path)):
            if not hasattr(filename, "name"):
                raise TypeError(f"Could not coerce filename from {filename}")
            else:
                filename = filename.name
        exists = os.path.exists(filename + cls.extension)
        inst = cls(filename)

        # File was empty
        if len(inst) == 0 and exists:
            exists = False
        return inst, exists

    @classmethod
    def exists(cls, filename: Union[str, pathlib.Path, io.FileIO]):
        if not isinstance(filename, (str, pathlib.Path)):
            if not hasattr(filename, "name"):
                raise TypeError(f"Could not coerce filename from {filename}")
            else:
                filename = filename.name
        exists = os.path.exists(filename + cls.extension)
        return exists

    def __init__(self, filename):
        self.filename = filename
        self.index_filename = self.filename + self.extension
        self._cache = None
        self.connect()
        self._size = len(self)
        self._size_uncommitted = 0

    def connect(self, create=None):
        filename = self.index_filename
        if os.path.exists(filename):
            if create:
                logger.debug(f'INFO: Deleting previous index file {filename}')
                os.remove(filename)
        logger.debug(f'INFO: Creating index file {filename}')
        engine = create_engine("sqlite:///"+filename)
        Base.metadata.create_all(engine)

        session = scoped_session(sessionmaker(bind=engine))
        self.session = session
        self.engine = engine
        self._cache = None

    def add(self, number, offset, name, analyte=None, attributes=None):
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
        record = SpectrumLibraryIndexRecord(number=number, offset=offset, name=name,
                                            index=self._size + self._size_uncommitted, analyte=analyte)
        self._size_uncommitted += 1
        if attributes is not None:
            raise NotImplementedError("Record attribute storage is not implemented")
        self.session.add(record)

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
        record = ClusterSpectrumLibraryIndexRecord(number=number, offset=offset)
        if attributes is not None:
            raise NotImplementedError("Record attribute storage is not implemented")
        self.session.add(record)

    def commit(self):
        """Persist any new entries to disk."""
        self._size += self._size_uncommitted
        self._size_uncommitted = 0
        self.session.commit()

    def iter_clusters(self) -> Iterator[IndexRecordBase]:
        """Iterate over cluster entries in the index."""
        for record in self.session.query(ClusterSpectrumLibraryIndexRecord).order_by(
                ClusterSpectrumLibraryIndexRecord.number).yield_per(10000):
            yield record

    def iter_spectra(self):
        """Iterate over spectrum entries in the index."""
        for record in self.session.query(SpectrumLibraryIndexRecord).order_by(
                SpectrumLibraryIndexRecord.number).yield_per(10000):
            yield record

    def __getitem__(self, i):
        return self._get_by_index(i)

    def _get_by_index(self, i: Union[int, slice]) -> Union[SpectrumLibraryIndexRecord, List[SpectrumLibraryIndexRecord]]:
        if isinstance(i, slice):
            records = self.session.query(SpectrumLibraryIndexRecord).slice(i.start, i.stop).all()
            if i.step:
                raise NotImplementedError()
            return records
        else:
            record = self.session.query(SpectrumLibraryIndexRecord).offset(i).limit(1).first()
            return record

    def __len__(self):
        value = self.session.query(func.count(SpectrumLibraryIndexRecord.id)).scalar()
        return value

    def search(self, i, **kwargs):
        if i is None and kwargs:
            # Executing attribute query
            raise NotImplementedError()
        if isinstance(i, numbers.Integral):
            if i < 0:
                i = len(self) + i
            if self._cache is not None and self._cache.number == i:
                return self._cache
            records = self.session.query(SpectrumLibraryIndexRecord).filter(SpectrumLibraryIndexRecord.number == i).all()

            if len(records) == 1:
                return records[0]
            elif len(records) == 0:
                raise IndexError(i)
            else:
                raise ValueError(f"Too many records found for spectrum number {i}")
        elif isinstance(i, slice):
            start = i.start or 0
            end = i.stop or float('inf')
            records = self.session.query(SpectrumLibraryIndexRecord).filter(
                SpectrumLibraryIndexRecord.number >= start,
                SpectrumLibraryIndexRecord.number < end).all()
            return records
        else:
            records = self.session.query(SpectrumLibraryIndexRecord).filter(
                SpectrumLibraryIndexRecord.name == i).all()
            if not records:
                raise KeyError(i)
            elif len(records) == 1:
                return records[0]
            else:
                return records

    def search_clusters(self, i, **kwargs):
        if i is None and kwargs:
            # Executing attribute query
            raise NotImplementedError()
        if isinstance(i, numbers.Integral):
            if i < 0:
                i = len(self) + i
            if self._cache is not None and self._cache.number == i:
                return self._cache
            records = self.session.query(ClusterSpectrumLibraryIndexRecord).filter(
                ClusterSpectrumLibraryIndexRecord.number == i).all()

            if len(records) == 1:
                return records[0]
            elif len(records) == 0:
                raise IndexError(i)
            else:
                raise ValueError(f"Too many records found for spectrum number {i}")
        elif isinstance(i, slice):
            start = i.start or 0
            end = i.stop or float('inf')
            records = self.session.query(ClusterSpectrumLibraryIndexRecord).filter(
                ClusterSpectrumLibraryIndexRecord.number >= start,
                ClusterSpectrumLibraryIndexRecord.number < end).all()
            return records
        else:
            raise NotImplementedError()
