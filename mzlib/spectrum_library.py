#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re
import timeit
import os

from mzlib.spectrum_library_index import SpectrumLibraryIndex
from mzlib.spectrum import Spectrum


debug = False

class SpectrumLibrary:
    """
    SpectrumLibrary - Class for a spectrum library

    Attributes
    ----------
    format : string
        Name of the format for the current encoding of the library.

    Methods
    -------
    read_header - Read just the header of the whole library
    read - Read the entire library into memory
    write - Write the library to disk
    create_index - Create an index file for this library
    transform - Not quite sure what this is supposed to be
    get_spectrum - Extract a single spectrum by identifier
    find_spectra - Return a list of spectra given query constraints

    """


    def __init__(self, identifier=None, name=None, filename=None, format=None):
        """
        __init__ - SpectrumLibrary constructor

        Parameters
        ----------
        format : string
            Name of the format for the current encoding of the library.

        """

        self.identifier = identifier
        self.name = name
        self.filename = filename
        self.format = format

        #### If we already have a filename, look for or create an index
        self.index = None
        if self.filename is not None:
            self.index = SpectrumLibraryIndex( library_filename=self.filename )


    #### Define getter/setter for attribute identifier
    @property
    def identifier(self):
        return(self._identifier)
    @identifier.setter
    def identifier(self, identifier):
        self._identifier = identifier

    #### Define getter/setter for attribute name
    @property
    def name(self):
        return(self._name)
    @name.setter
    def name(self, name):
        self._name = name

    #### Define getter/setter for attribute filename
    @property
    def filename(self):
        return(self._filename)
    @filename.setter
    def filename(self, filename):
        self._filename = filename

    #### Define getter/setter for attribute format
    @property
    def format(self):
        return(self._format)
    @format.setter
    def format(self, format):
        self._format = format



    def read_header(self):
        """
        read_header - Read just the header of the whole library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        filename = self.filename
        if debug: eprint(f"INFO: Reading library header from {filename}")
        if filename is None:
            eprint("ERROR: Unable to read library with no filename")
            return(False)
        with open(filename, 'r') as stream:
            first_line = stream.readline()
            if re.match("Name: ",first_line):
                if debug: eprint("INFO: This appears to be a headerless MSP file")
                self.format = "msp"
                self.header = []
                return(True)
        return(False)


    def read(self, create_index=None):
        """
        read - Read the entire library into memory

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Check that the spectrum library filename isvalid
        filename = self.filename
        if debug: eprint(f'INFO: Reading spectra from {filename}')
        if filename is None:
            eprint("ERROR: Unable to read library with no filename")
            return(False)

        #### If an index hasn't been opened, open it now
        if create_index is not None:
            if self.index is None:
                self.index = SpectrumLibraryIndex( library_filename=self.filename )
                self.index.create_index()

        #### Determine the filesize
        file_size = os.path.getsize(filename)
        if debug: eprint(f"INFO: File size is {file_size}")

        with open(filename, 'r') as infile:
            state = 'header'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''
            if debug: eprint("INFO: Reading..",end='',flush=True)
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                #file_offset = infile.tell()
                file_offset += len(line)

                line = line.rstrip()
                if state == 'header':
                    if re.match('Name: ',line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('Name: ',line):
                        if len(spectrum_buffer) > 0:
                            #parse(spectrum_buffer)
                            if create_index is not None:
                                self.index.add_spectrum( number=n_spectra + start_index, offset=spectrum_file_offset, name=spectrum_name, peptide_sequence=None )
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            if int(n_spectra/1000) == n_spectra/1000:
                                self.index.commit()
                                percent_done = int(file_offset/file_size*100+0.5)
                                eprint(str(percent_done)+"%..",end='',flush=True)

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match('Name:\s+(.+)',line).group(1)
                        #print(spectrum_name)
                    spectrum_buffer.append(line)
                #if n_spectra > 5:
                #    break

            #### Process the last spectrum in the buffer
            #parse(spectrum_buffer)
            if create_index is not None:
                self.index.add_spectrum( number=n_spectra + start_index, offset=spectrum_file_offset, name=spectrum_name, peptide_sequence=None )
                self.index.commit()
            n_spectra += 1
            if debug:
                eprint()
                eprint(f"INFO: Read {n_spectra} spectra from {filename}")

            #### Flush the index
            #self.index.commit()
        return(n_spectra)


    def read_spectrum(self, offset=None):
        """
        read - Read the entire library into memory

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Check that an offset is supplied
        if offset is None:
            eprint("ERROR: Required parameter offset is not supplied")
            return(False)

        #### Check that the spectrum library filename is valid
        filename = self.filename
        if debug: eprint(f'INFO: Reading spectrum from {filename} at offset {offset}')
        if filename is None:
            eprint("ERROR: Unable to read library with no filename")
            return(False)

        with open(filename, 'r') as infile:
            infile.seek(offset)
            state = 'body'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''
            for line in infile:
                line_beginning_file_offset = file_offset
                file_offset += len(line)
                line = line.rstrip()
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('Name: ',line):
                        if len(spectrum_buffer) > 0:
                            #parse(spectrum_buffer)
                            return(spectrum_buffer)
                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match('Name:\s+(.+)',line).group(1)
                        #print(spectrum_name)
                    spectrum_buffer.append(line)

            #### We will end up here if this is the last spectrum in the file
            #parse(spectrum_buffer)
            return(spectrum_buffer)

        return(False)


    def write(self):
        """
        write - Write the library to disk

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
        create_index - Create an index file for this library

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        self.read(create_index=True)
        return()


    def get_spectrum(self,spectrum_index_number=None,spectrum_name=None):
        """
        get_spectrum - Extract a single spectrum by identifier

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here
        if self.index is None:
            self.index = SpectrumLibraryIndex( library_filename=self.filename )

        #### If spectrum_index_number was specified, find the spectrum by that
        if spectrum_index_number is not None:
            offset = self.index.get_offset(spectrum_index_number=spectrum_index_number)
            if offset is not None:
                if debug: print(f'Found offset {offset} for spectrum {spectrum_index_number}')
            else:
                if debug: print(f'Unable to find offset for spectrum {spectrum_index_number}')
            spectrum_buffer = self.read_spectrum(offset=offset)
            return(spectrum_buffer)
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

    #### Create a new RTXFeedback object
    spectrum_library = SpectrumLibrary()
    spectrum_library.filename = "../refData/sigmaups1_consensus_final_true_lib.msp"
    spectrum_library.read_header()
 
    spectrum_buffer = spectrum_library.get_spectrum(spectrum_index_number=2000)
    spectrum = Spectrum()
    spectrum.parse(spectrum_buffer)
    buffer = spectrum.write(format="text")
    print(buffer)
    print()

    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example()
    return()


if __name__ == "__main__": main()

