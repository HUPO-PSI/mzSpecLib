#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import os
import os.path
import argparse
import re

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../lib")
from SpectrumLibraryCollection import SpectrumLibraryCollection

def main():

    argparser = argparse.ArgumentParser(description='Front-end for maintenance of a spectral library collection')

    argparser.add_argument('--list', action='store_true', help="List the available libraries in the collection")
    argparser.add_argument('--update', action='store_true', help="Rescan the spectralLibraries directory and update the collection")
    argparser.add_argument('--DROP', action='store_true', help="DROP the current collection database and re-create it. Careful!")

    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Ensure that library_file was passed
    if params.list is False and params.update is False:
        print("ERROR: Nothing to do. See --help for more information")
        return()

    collection_dir = os.path.dirname(os.path.abspath(__file__))+"/../spectralLibraries"

    #### If --list, then list existing libraries
    if params.DROP is True:
        if os.path.exists(collection_dir + "/SpectrumLibraryCollection.sqlite"):
            print("INFO: DELETing the collections database and re-creating")
            os.remove(collection_dir + "/SpectrumLibraryCollection.sqlite")

    #### Get a listing of what entries are already in the collection
    spectrum_library_collection = SpectrumLibraryCollection(collection_dir + "/SpectrumLibraryCollection.sqlite")
    libraries = spectrum_library_collection.get_libraries()

    #### If --update, then sync the files in spectralLibraries with the collection
    if params.update is True:

        #### Also create a dict for the filenames
        libraries_dict = {}
        counter = 0
        for library in libraries:
            libraries_dict[library.original_name] = counter
            counter += 1

        #### Read the metadata file for the collection
        metadatafile = collection_dir + "/SpectrumLibraryCollection.tsv"
        metadata_dict = {}
        metadata_files = []
        with open(metadatafile, 'r') as infile:
            print(f"INFO: Reading {metadatafile}")
            for line in infile:
                line = line.rstrip()
                columns = line.split("\t")
                if columns[4] == "filename":
                    continue
                metadata_dict[columns[4]] = columns
                metadata_files.append(columns[4])

        #### Loop over all files in the directory and put them in a hash
        local_file_dict = {}
        for filename in os.listdir(collection_dir):
            match = re.fullmatch("(.+)\.(msp|sptxt)",filename)
            if match is not None:
                local_file_dict[filename] = 1

        #### Loop over all files in the metadata files and check against the database
        for filename in metadata_files:
            if filename in local_file_dict:
                print(f"Processing {filename}")

                #### Check to see if there is an index file yet
                if os.path.exists(f"{collection_dir}/{filename}.splindex"):
                    print("  Index already exists")
                else:
                    print("  Need to create an index")
                    os.system(f"{collection_dir}/../bin/index_library.py --lib {filename}")

                #### Check to see if the file is already registered
                if filename in libraries_dict:
                    print(f"  This library is already in the collection as {libraries[libraries_dict[filename]].id_name}")
                    print(f"    (stored version={libraries[libraries_dict[filename]].version}")
                    version = metadata_dict[filename][3].strip('"')
                    print(f"    (metadata sheet version={version}")
                    spectrum_library_collection.update_library_metadata(id=libraries[libraries_dict[filename]].library_record_id, version=version)

                else:
                    print(f"  Need to create a record for filename")
                    if filename in metadata_dict:
                        version = metadata_dict[filename][3]
                        spectrum_library_collection.add_library(original_name=filename, version=version)
                    else:
                        print(f"ERROR: There is no metadata entry for {filename} yet. Please create one and update again")
                        return()

            else:
                print(f"ERROR: filename '{filename}' in the metadata file not found locally")
                return()

    #### If --list, then list existing libraries
    if params.list is True:
        libraries = spectrum_library_collection.get_libraries()
        if ( len(libraries) == 0 ):
            print("The library collection is empty")
        else:
            for library in libraries:
                print("\t".join([str(library.library_record_id),library.id_name,str(library.version),library.original_name]))


    return()

if __name__ == "__main__": main()
