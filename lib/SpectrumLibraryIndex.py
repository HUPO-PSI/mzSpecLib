#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import os
from datetime import datetime

from sqlalchemy import Column, ForeignKey, Integer, Float, String, DateTime, Text, PickleType, LargeBinary
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import desc
from sqlalchemy import inspect

Base = declarative_base()

debug = True


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
  peptide_sequence = Column(String(2014), nullable=True)


class SpectrumLibraryIndex:
    """
    SpectrumLibraryIndex - Class for a spectrum library index

    Attributes
    ----------
    columns : array
        Names of the columns in the data matrix

    Methods
    -------
    get_offset - Get the offset for a spectrum in the library based on the spectrum_index or spectrum_name
    find_offsets - Return an array of offsets of spectra that match the input parameters
    create_index - Create a new index for a library
    add_spectrum - Add a spectrum to the index

    """


    #### Constructor
    def __init__(self, library_filename=None, version=None, n_spectra=None, library_datetime=None, columns=None):
        """
        __init__ - SpectrumLibraryIndex constructor

        Parameters
        ----------
        columns : array
            Names of the columns in the data matrix

        """

        self.library_filename = library_filename
        self.version = "0.1"
        self.n_spectra = 0
        self.library_datetime = None
        self.columns = [ 'number', 'offset', 'name', 'peptide_sequence' ]
        self.status = 'closed'
        self.uncommitted_transactions = 0

        if library_filename is None:
            raise Exception('Library filename missing')
        filename = self.library_filename + '.splindex'
        if filename is None:
            raise Exception('Missing library_filename')
        if os.path.exists(filename):
            self.connect()
            self.status = 'OK'
        else:
            self.create_database()
            self.status = 'OK'


    #### Destructor
    def __del__(self):
        if self.library_filename is not None:
            self.disconnect()


    #### Define getter/setter for attribute library_filename
    @property
    def library_filename(self):
        return(self._library_filename)
    @library_filename.setter
    def library_filename(self, library_filename):
        self._library_filename = library_filename

    #### Define getter/setter for attribute version
    @property
    def version(self):
        return(self._version)
    @version.setter
    def version(self, version):
        self._version = version

    #### Define getter/setter for attribute n_spectra
    @property
    def n_spectra(self):
        return(self._n_spectra)
    @n_spectra.setter
    def n_spectra(self, n_spectra):
        self._n_spectra = n_spectra

    #### Define getter/setter for attribute library_datetime
    @property
    def library_datetime(self):
        return(self._library_datetime)
    @library_datetime.setter
    def library_datetime(self, library_datetime):
        self._library_datetime = library_datetime

    #### Define getter/setter for attribute columns
    @property
    def columns(self):
        return(self._columns)
    @columns.setter
    def columns(self, columns):
        self._columns = columns


    #### Define attribute session
    @property
    def session(self) -> str:
        return self._session

    @session.setter
    def session(self, session: str):
        self._session = session


    #### Define attribute engine
    @property
    def engine(self) -> str:
        return self._engine

    @engine.setter
    def engine(self, engine: str):
        self._engine = engine


    #### Delete and create the database. Careful!
    def create_database(self):
        filename = self.library_filename + '.splindex'
        if os.path.exists(filename):
            if debug: eprint(f'INFO: Deleting previous index file {filename}')
            os.remove(filename)
        if debug: eprint(f'INFO: Creating index file {filename}')
        engine = create_engine("sqlite:///"+filename)
        Base.metadata.create_all(engine)

        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        self.session = session
        self.engine = engine

        index_attribute = SpectrumLibraryIndexAttribute( name='version', value=self.version )
        session.add(index_attribute)
        index_attribute = SpectrumLibraryIndexAttribute( name='n_spectra', value=0 )
        session.add(index_attribute)

        session.flush()
        session.commit()


    #### Create and store a database connection
    def connect(self):
        filename = self.library_filename + '.splindex'
        if debug: eprint(f'INFO: Opening index file {filename}')
        engine = create_engine("sqlite:///"+filename)
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        self.session = session
        self.engine = engine


    #### Destroy the database connection
    def disconnect(self):
        filename = self.library_filename + '.splindex'
        if self.uncommitted_transactions > 0:
            if debug: eprint(f'INFO: Committing open transactions')
            self.session.flush()
            self.session.commit()
            self.uncommitted_transactions = 0

        if debug: eprint(f'INFO: Closing index file {filename}')
        session = self.session
        engine = self.engine
        session.close()
        engine.dispose()
        self.session = None
        self.engine = None


    #### Destroy the database connection
    def commit(self):
        if self.uncommitted_transactions > 0:
            #if debug: eprint(f'INFO: Committing open transactions')
            self.session.flush()
            self.session.commit()
            self.uncommitted_transactions = 0


    def get_offset(self, spectrum_index_number=None, spectrum_name=None):
        """
        get_offset - Get the offset for a spectrum in the library based on the spectrum_index or spectrum_name

        Extended description of function.

        Parameters
        ----------
        spectrum_index : integer
            Index number of the spectrum to select
        spectrum_name : string
            Name of the spectrum to select

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if spectrum_index_number is not None:
            try:
                records = self.session.query(SpectrumLibraryIndexRecord).filter(SpectrumLibraryIndexRecord.number==spectrum_index_number).all()
            except:
                session.rollback()
                raise
            if len(records) > 1:
                raise Exception('Too many records')
            return(records[0].offset)

        return()


    def find_offsets(self):
        """
        find_offsets - Return an array of offsets of spectra that match the input parameters

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here

        return()


    def create_index(self):
        """
        create_index - Create a new index for a library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if self.session is not None:
            self.disconnect()
        self.create_database()
        return(True)


    def add_spectrum(self, number=None, offset=None, name=None, peptide_sequence=None):
        """
        add_spectrum - Add a spectrum to the index

        Extended description of function.

        Parameters
        ----------
        number : integer
            Index number of the spectrum to add
        offset : integer
            File offset of the spectrum to add
        name : string
            Name of the spectrum to add
        peptide_sequence : string
            Unmodified peptide sequence of the spectrum to add

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        #if debug: eprint("INFO: Adding an index entry")
        session = self.session
        index_record = SpectrumLibraryIndexRecord( number=number, offset=offset, name=name, peptide_sequence=peptide_sequence )
        session.add(index_record)
        self.uncommitted_transactions += 1
        if self.uncommitted_transactions >= 5000:
            session.flush()
            session.commit()
            self.uncommitted_transactions = 0

        self.status = 'uncommitted changes'
        return()



#### Example using this class
def example():

    #### Create a new RTXFeedback object
    index = SpectrumLibraryIndex(library_filename='../refData/sigmaups1_consensus_final_true_lib.msp')
    print(index.version)
    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example()
    return()


if __name__ == "__main__": main()

