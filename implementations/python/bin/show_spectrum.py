#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import os
import argparse
import os.path
import timeit

basedir = "some hard coded location"
basedir = os.path.dirname(os.path.abspath(__file__))+"/.."

sys.path.append(basedir+"/lib")
from mzlib import SpectrumLibraryCollection
from mzlib import SpectrumLibrary
from mzlib import Spectrum
from mzlib import UniversalSpectrumIdentifier


def main():

    argparser = argparse.ArgumentParser(description='Reads one spectrum from a file and prints to stdout')

    argparser.add_argument('--library_file', action='store', help="Name of the library file to access")
    argparser.add_argument('--index_number', action='store', help="Index number of the spectrum to display")
    argparser.add_argument('--usi', action='store', help="Universal Spectrum Identifier of the spectrum to display")
    argparser.add_argument('--output_format', action='store', default='text', help="Format use when writing the spectrum (one of 'text', 'json', 'msp')")

    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #print("Content-type: text/plain\n")
    #print(os.environ)

    library_file = params.library_file
    index_number = params.index_number

    #### Ensure that either a USI or a library_file and index_number was passed
    if params.usi is None or params.usi == "":
        if params.library_file is None or params.library_file == "":
            print("ERROR: Parameter --usi or --library_file must be provided. See --help for more information")
            return()

        if params.index_number is None or params.index_number == "":
            print("ERROR: Parameter --usi or index_number must be provided. See --help for more information")
            return()

    #### If there was a USI, then parse it
    else:
        usi = UniversalSpectrumIdentifier(params.usi)
        if not usi.is_valid:
            print(f"ERROR: {usi.error_code}: {usi.error_message}")
            return()
        if usi.dataset_identifier.startswith("PXL"):
            spec_lib_collection = SpectrumLibraryCollection(basedir + "/spectralLibraries/SpectrumLibraryCollection.sqlite")
            #print(f"Looking up library for {usi.dataset_identifier}")
            try:
                library = spec_lib_collection.get_library(identifier=usi.dataset_identifier, version=usi.ms_run_name)
            except Exception as error:
                print("ERROR:",error)
                return()
            #print("Found record: " + "\t".join([str(library.library_record_id),library.id_name,library.version,library.original_name]))
            library_file = basedir + "/spectralLibraries/" + library.original_name
            index_number = usi.index

    if not os.path.isfile(library_file):
        eprint(f"ERROR: File '{library_file}' not found or not a file")
        return()

    os.remove(library_file + '.splindex')

    spectrum_library = SpectrumLibrary(filename=library_file)
    spectrum_library.read(create_index=True)

    spectrum_buffer = spectrum_library.get_spectrum(spectrum_index_number=index_number)
    spectrum = Spectrum()
    spectrum.parse(spectrum_buffer, spectrum_index=index_number)
    buffer = spectrum.write(format=params.output_format)
    print(buffer)

    return()

if __name__ == "__main__": main()
