import os
import numbers
import pathlib
import logging

from sqlalchemy import Column, ForeignKey, Integer, Float, String, DateTime, Text, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker

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
    number = Column(Integer, nullable=False)
    offset = Column(Integer, nullable=False)
    name = Column(String(1024), nullable=False)
    analyte = Column(String(2014), nullable=True)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.number}, {self.offset}, {self.name}, {self.analyte})"


class SQLIndex(IndexBase):
    extension = '.splindex'

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
    def exists(cls, filename):
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
        self.connect()

    def connect(self, create=None):
        filename = self.index_filename
        if os.path.exists(filename):
            if create:
                logger.debug(f'INFO: Deleting previous index file {filename}')
                os.remove(filename)
        logger.debug(f'INFO: Creating index file {filename}')
        engine = create_engine("sqlite:///"+filename)
        Base.metadata.create_all(engine)

        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        self.session = session
        self.engine = engine

    def add(self, number, offset, name, analyte, attributes=None):
        record = SpectrumLibraryIndexRecord(number=number, offset=offset, name=name, analyte=analyte)
        if attributes is not None:
            raise NotImplementedError("Record attribute storage is not implemented")
        self.session.add(record)

    def commit(self):
        self.session.commit()

    def __iter__(self):
        for record in self.session.query(SpectrumLibraryIndexRecord).order_by(SpectrumLibraryIndexRecord.number):
            yield record

    def __getitem__(self, i):
        return self.search(i)

    def __len__(self):
        value = self.session.query(func.count(SpectrumLibraryIndexRecord.id)).scalar()
        return value

    def search(self, i, **kwargs):
        if i is None and kwargs:
            # Executing attribute query
            raise NotImplementedError()
        if isinstance(i, numbers.Integral):
            records = self.session.query(SpectrumLibraryIndexRecord).filter(SpectrumLibraryIndexRecord.number == i).all()
            if len(records) == 1:
                return records[0]
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

