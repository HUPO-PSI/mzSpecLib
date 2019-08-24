#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import os
import argparse
import os.path
import timeit

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from SpectrumLibrary import SpectrumLibrary
from Spectrum import Spectrum

def main():

    argparser = argparse.ArgumentParser(description='Reads one spectrum from a file and prints to stdout')

    argparser.add_argument('--library_file', action='store', help="Name of the library file to access")
    argparser.add_argument('--index_number', action='store', help="Index number of the spectrum to display")
    argparser.add_argument('--output_format', action='store', default='text', help="Format use when writing the spectrum (one of 'text', 'json', 'msp')")

    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Ensure that library_file was passed
    if params.library_file is None or params.library_file == "":
        print("ERROR: Parameter --library_file must be provided. See --help for more information")
        return()

    #### Ensure that library_file was passed
    if params.index_number is None or params.index_number == "":
        print("ERROR: Parameter --index_number must be provided. See --help for more information")
        return()

    if not os.path.isfile(params.library_file):
        eprint(f"ERROR: File '{params.library_file}' not found or not a file")
        return()

    spectrum_library = SpectrumLibrary()
    spectrum_library.filename = params.library_file

    spectrum_buffer = spectrum_library.get_spectrum(spectrum_index_number=params.index_number)
    spectrum = Spectrum()
    spectrum.parse(spectrum_buffer)
    buffer = spectrum.write(format=params.output_format)
    print(buffer)
    print()

    return()

if __name__ == "__main__": main()
