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

debug = False

#### Define the database tables as classes
class LibraryRecord(Base):
  __tablename__ = 'library_record'
  library_record_id = Column(Integer, primary_key=True)
  id_name = Column(String(255), nullable=False)
  version = Column(String(255), nullable=False)
  status = Column(String(255), nullable=False)
  original_name = Column(String(255), nullable=False)
  original_checksum = Column(String(255), nullable=True)
  n_entries = Column(Integer, nullable=True)
  record_datetime = Column(DateTime, nullable=False)


class SpectrumLibraryCollection:
    """
    SpectrumLibraryCollection - Class for a collection of spectrum libraries

    Attributes
    ----------
    filename : string
        Filename of the SQLite database file that contains information about the collection of libraries available.

    Methods
    -------
    create - Create a new spectral library collection
    show - Return a string that summarizes the state of the collection
    get_libraries - Return a list of available libraries
    get_library - Return attributes of a specific library
    add_library - Add a new library
    create_index - Create a master index from all the constituent library indexes to be able to find spectra in any library
    find_spectra - Return a list of spectra given query constraints

    """


    def __init__(self, filename=None):
        """
        __init__ - SpectrumLibraryCollection constructor

        Parameters
        ----------
        filename : string
            Filename of the SQLite database file that contains information about the collection of libraries available.

        """

        self.filename = filename
        if os.path.exists(self.filename):
            if debug: eprint(f"INFO: Library collection {self.filename} exists")
            self.connect()
        else:
            if debug: eprint(f"INFO: Library collection {self.filename} not found. Will create new.")
            self.createDatabase()


    #### Destructor
    def __del__(self):
        self.disconnect()


    #### Define getter/setter for attribute filename
    @property
    def filename(self):
        return(self._filename)
    @filename.setter
    def filename(self, filename):
        self._filename = filename


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
    def createDatabase(self):
        if os.path.exists(self.filename):
            eprint("INFO: Deleting previous database file " + self.filename)
            os.remove(self.filename)
        if debug: eprint("INFO: Creating database " + self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        Base.metadata.create_all(engine)
        self.connect()


    #### Create and store a database connection
    def connect(self):
        if debug: eprint("INFO: Connecting to database " + self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        self.session = session
        self.engine = engine


    #### Destroy the database connection
    def disconnect(self):
        if debug: eprint("INFO: Disconnecting from database " + self.filename)
        session = self.session
        engine = self.engine
        session.close()
        engine.dispose()


    def create(self, overwrite_existing=False):
        """
        create - Create a new spectral library collection

        Extended description of function.

        Parameters
        ----------
        overwrite_existing : boolean
            Set to true in order to write over the previous file if it exists

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug: eprint("INFO: Creating database " + self.filename)
        if os.path.exists(self.filename):
            os.remove(self.filename)
        engine = create_engine("sqlite:///"+self.filename)
        Base.metadata.create_all(engine)
        self.connect()
        return()


    def get_libraries(self):
        """
        get_libraries - Return a list of available libraries

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug: eprint("INFO: Fetching all libraries")
        session = self.session
        libraries = session.query(LibraryRecord).all()
        return(libraries)


    def get_library(self,identifier=None,version=None,filename=None):
        """
        get_library - Return attributes of a specific library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        session = self.session
        if identifier is not None and identifier > "":
            libraries = session.query(LibraryRecord).filter(LibraryRecord.id_name==identifier).order_by(desc(LibraryRecord.version)).all()
            if len(libraries) == 0:
                raise Exception(f"No library with identifier {identifier} found")
            else:
                libraryListStr = ""
                for library in libraries:
                    libraryListStr += library.version+","
                    if version is not None and version > "":
                        if version == library.version:
                            return(library)
                if version is not None and version > "":
                    raise Exception(f"Unable to find version {version} of library {identifier}")
            if len(libraries) == 1:
                return(library)
            raise Exception(f"There are several version of this library ({libraryListStr}). Please specify a version")
            return()

        elif filename is not None and filename > "":
            raise Exception("Search by filename not implemented")
        else:
            raise Exception("Not enough information to find library")

        return()


    def add_library(self, original_name, version="1"):
        """
        add_library - Add a new library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if debug: eprint("INFO: Adding a library entry")
        session = self.session
        library_record = LibraryRecord(id_name="PXL000000", version=version,
            status="initial_add",original_name=original_name,record_datetime=datetime.now())
        session.add(library_record)
        session.flush()
        assert(library_record.library_record_id)
        idstr = str(library_record.library_record_id)
        if debug: eprint(f"INFO: Returned id={idstr}")
        idstr_length = len(idstr)
        assert(idstr_length)
        padding = "000000"
        new_idstr = "PXL" + padding[0:len(padding)-idstr_length] + idstr
        library_record.id_name = new_idstr
        session.flush()
        session.commit()
        return()


    def create_index(self):
        """
        create_index - Create a master index from all the constituent library indexes to be able to find spectra in any library

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


    def find_spectra(self):
        """
        find_spectra - Return a list of spectra given query constraints

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


#### Example using this class
def example():

    #### Create a new or attach to existing library collection
    spec_lib_collection = SpectrumLibraryCollection("zzTest.sqlite")

    #### Add a library
    spec_lib_collection.add_library()

    #### Show all libraries
    libraries = spec_lib_collection.get_libraries()
    if ( libraries is None ):
        print("The library collection is empty")
    else:
        for library in libraries:
            print(library.library_record_id,library.id_name,library.original_name)

    return()


#### Example using this class
def example2():

    #### Create a new or attach to existing library collection
    spec_lib_collection = SpectrumLibraryCollection("../spectralLibraries/SpectrumLibraryCollection.sqlite")

    try:
        library = spec_lib_collection.get_library(identifier="PXL000003", version="05-24-2011")
    except Exception as error:
        print("ERROR:",error)
        return()
    print("\t".join([str(library.library_record_id),library.id_name,library.version,library.original_name]))
    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example2()
    return()


if __name__ == "__main__": main()

