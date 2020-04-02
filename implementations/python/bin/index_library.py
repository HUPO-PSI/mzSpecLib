#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import os
import argparse
import os.path
import timeit

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/..")
from mzlib import SpectrumLibrary

def main():

    argparser = argparse.ArgumentParser(description='Creates an index for an MSP spectral library file')

    argparser.add_argument('--library_file', action='store',
        help='Name of the library to index')

    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Ensure that library_file was passed
    if params.library_file is None or params.library_file == "":
        print("ERROR: Parameter --library_file must be provided. See --help for more information")
        return()

    if not os.path.isfile(params.library_file):
        eprint(f"ERROR: File '{params.library_file}' not found or not a file")
        return()

    spectrum_library = SpectrumLibrary()
    spectrum_library.filename = params.library_file

    t0 = timeit.default_timer()
    spectrum_library.create_index()
    t1 = timeit.default_timer()
    print()
    print('INFO: Elapsed time: ' + str(t1-t0))

if __name__ == "__main__": main()
