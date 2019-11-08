#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re
from ontology_term import OntologyTerm


#############################################################################
#### Ontology class
class Ontology(object):


    #########################################################################
    #### Constructor
    def __init__(self, filename=None, verbose=0):

        self.filename = filename
        self.verbose = verbose

        self.is_valid = False
        self.n_terms = 0
        self.primary_prefix = None
        self.prefixes = {}
        self.header_line_list = []
        self.other_line_list = []
        self.term_list = []
        self.terms = {}
        self.names = {}
        self.uc_names = {}
        self.n_errors = 0
        self.error_code = None
        self.error_message = None

        self.uc_search_string = None

        #### If we have been given a filename on construction, read it right away
        if filename:
            self.read()
        

    #########################################################################
    #### parse the file
    def read(self, filename=None, verbose=0):
        verboseprint = print if verbose>0 else lambda *a, **k: None

        #### Determine the filename to read
        if filename is not None:
            self.filename = filename
        filename = self.filename

        #### Set up some beginning statement
        state = 'header'
        terms_list = []
        terms = {}
        current_term = []

        verboseprint(f"INFO: Reading file '{filename}'")
        with open(filename, encoding="latin-1", errors="replace") as infile:
            for line in infile:
                line = line.rstrip()

                #### Process the header
                if state == 'header':
                    match = re.search("^\s*\[Term\]\s*$",line)
                    if match:
                        state = 'term'
                    else:
                        self.header_line_list.append(line)

                #### Process the other elements in the file
                if state == 'other':
                    match = re.search("^\s*\[Term\]\s*$",line)
                    if match:
                        state = 'term'
                    else:
                        self.other_line_list.append(line)

                #### Process the term section
                if state == 'term':

                    #### Skip an empty line
                    match = re.search("^\s*$",line)
                    if match:
                        continue

                    #### If this is a new element
                    match = re.search("^\s*\[",line)
                    if match:

                        #### If this is a new [Term]
                        match = re.search("^\s*\[Term\]\s*$",line)
                        if match:

                            #### If there is currently something in the buffer, process it
                            #### WARNING: If any changes are made here, you also need to update after the loop for the processing of the last term
                            if len(current_term) > 0:
                                #### Store this term
                                term = OntologyTerm(line_list=current_term, verbose=verbose)
                                if term.is_obsolete is False:
                                    self.term_list.append(term.curie)
                                    if term.curie in self.terms:
                                        set_error("Duplicate term!")
                                    else:
                                        self.terms[term.curie] = term
                                        if term.prefix not in self.prefixes:
                                            self.prefixes[term.prefix] = 0
                                        self.prefixes[term.prefix] += 1

                                #### Clear current term data and update counters
                                current_term = []
                                self.n_terms += 1

                            #### Append the current line to the working term buffer
                            current_term.append(line)

                        #### Otherwise this is the start of some non-Term thing
                        else:
                            state = 'other'
                            self.other_line_list.append(line)

                    #### Otherwise, just append the current line to the working term buffer
                    else:
                        current_term.append(line)

        #### Process a last term that still may be in the buffer
        if len(current_term) > 0:
            #### Store this term
            term = OntologyTerm(line_list=current_term, verbose=verbose)
            if term.is_obsolete is False:
                self.term_list.append(term.curie)
                if term.curie in self.terms:
                    set_error("Duplicate term!")
                else:
                    self.terms[term.curie] = term
                    if term.prefix not in self.prefixes:
                        self.prefixes[term.prefix] = 0
                    self.prefixes[term.prefix] += 1

            current_term = []
            self.n_terms += 1

        #### Now map the parentage structure into children
        self.map_children(verbose=verbose)
        #### And create the map of names
        self.create_name_map(verbose=verbose)

        #### Set the is_valid state
        if self.n_errors == 0:
           self.is_valid = True

        else:
            self.is_valid = False
            verboseprint(f"Number of errors: {self.n_error}")
 
 
    #########################################################################
    #### Map all the parent relationships to child relationships for the parent
    def map_children(self, verbose=0):
        verboseprint = print if verbose>0 else lambda *a, **k: None

        verboseprint(f"INFO: Mapping parents to children")
        for curie in self.term_list:
            term = self.terms[curie]
            parents = term.parents
            for parent in parents:
                parent_curie = parent['curie']
                type = parent['type']
                new_type = '??'
                if type == 'is_a': new_type = 'has_subclass'
                if type == 'part_of': new_type = 'has_part'
                if parent_curie in self.terms:
                    self.terms[parent_curie].children.append( { 'type': new_type, 'curie': curie } )
                else:
                    if parent_curie != 'UO:0000000':
                        verboseprint(f"{curie} has parent {parent_curie}, but this curie is not found in this ontology")


    #########################################################################
    #### Create a dict of all the names and synonyms
    def create_name_map(self, verbose=0):
        verboseprint = print if verbose>0 else lambda *a, **k: None

        verboseprint(f"INFO: Creating a dict of all names and synonyms")
        for curie in self.term_list:
            term = self.terms[curie]
            names = [ term.name ]
            for synonym in term.synonyms:
                names.append(synonym['term'])
            for name in names:
                if name is None:
                    #print(f"WARNING: Term {curie} has no name!")
                    name = curie
                if name in self.names:
                    self.names[name].append(curie)
                else:
                    self.names[name] = [ curie ]
                uc_name = name.upper()
                if uc_name in self.uc_names:
                    self.uc_names[uc_name].append(curie)
                else:
                    self.uc_names[uc_name] = [ curie ]


    #########################################################################
    #### Get a list of all children of a term
    def get_children(self, parent_curie, return_type='ucdict'):

        if parent_curie not in self.terms:
            return([])

        children_curies = {}

        parent_term = self.terms[parent_curie]
        children = parent_term.children

        while len(children) > 0:
            new_children = []
            for child in children:
                child_curie = child['curie']
                children_curies[child_curie] = 1
                child_term = self.terms[child_curie]
                if len(child_term.children) > 0:
                    new_children.extend(child_term.children)
            children = []
            children.extend(new_children)

        result_list = []
        results = {}
        for child in children_curies:
            if return_type == 'ucdict':
                results[self.terms[child].name.upper()] = [child]
            if return_type == 'uclist':
                result_list.append(self.terms[child].name.upper())

        if return_type == 'ucdict':
            return(results)
        if return_type == 'uclist':
            return(result_list)
        


    #########################################################################
    #### Fuzzy search for a string
    def fuzzy_search(self, search_string, max_hits=15, children_of=None):

        match_term_list = []
        match_curies = {}

        eprint(f"INFO: Executing fuzzy search for '{search_string}'")
        search_space = self.uc_names
        if children_of is not None:
            search_space = self.get_children(parent_curie=children_of, return_type='ucdict')

        self.uc_search_string = search_string.upper()
        match_list = filter(self.filter_starts_with,search_space)
        for match in match_list:
            curies = search_space[match]
            curie = curies[0]
            if curie in match_curies: continue
            match_curies[curie] = 1
            #print("--",curie)
            term = { 'curie': curie, 'name': self.terms[curie].name, 'sort': 1 }
            match_term_list.append(term)

        count = len(match_term_list)

        if count < max_hits:
            matches = filter(self.filter_contains,search_space)
            for match in matches:
                curies = search_space[match]
                curie = curies[0]
                if curie in match_curies: continue
                match_curies[curie] = 1
                #print("==",curie)
                term = { 'curie': curie, 'name': self.terms[curie].name, 'sort': 2 }
                match_term_list.append(term)

        sorted_match_term_list = sorted(match_term_list,key=sort_by_relevance)
        if len(sorted_match_term_list) > max_hits:
            del sorted_match_term_list[max_hits:]

        for match in sorted_match_term_list:
            del match['sort']

        return(sorted_match_term_list)


    #########################################################################
    # Print out some information about the ontology
    def show(self):
        print(f"filename: {self.filename}")
        print(f"is_valid: {self.is_valid}")
        print(f"Number of errors: {self.n_errors}")
        print(f"Number of terms: {self.n_terms}")
        print(f"Number of header lines: {len(self.header_line_list)}")
        print(f"prefixes: {self.prefixes}")

        print(f"Number of other lines: {len(self.other_line_list)}")
        if len(self.other_line_list) > 0:
            i_other_line = 0
            for line in self.other_line_list:
                print("  >"+line)
                i_other_line += 1
                if i_other_line > 15: break
            print("  >...")


    #########################################################################
    #### filtering routines
    def filter_starts_with(self,x):
        match = re.match(self.uc_search_string,x)
        if match: return(True)
        return(False)
    def filter_contains(self,x):
        match = re.search(self.uc_search_string,x)
        if match: return(True)
        return(False)


#########################################################################
#### sorting routines
def sort_by_relevance(x):
    value = x['sort'] * 1000 + len(x['name'])
    return(value)


#########################################################################
#### A very simple example of using this class
def psims_example():
    ontology = Ontology(filename='psi-ms.obo',verbose=1)
    ontology.show()
    print("============================")
    term = ontology.terms["MS:1002286"]
    term.show()
    print("============================")
    name = 'QIT'
    print(f"Searching for '{name}'")
    if name in ontology.names:
        curies = ontology.names[name]
        for curie in curies:
            term = ontology.terms[curie]
            term.show()
    print("============================")
    name = 'bit'
    result_list = ontology.fuzzy_search(search_string=name)
    for item in result_list:
        print(item)
    print("============================")
    name = 'bit'
    result_list = ontology.fuzzy_search(search_string=name,children_of="MS:1000031")
    for item in result_list:
        print(item)


#########################################################################
#### A very simple example of using this class
def po_example():
    ontology = Ontology(filename='plant-ontology.obo',verbose=1)
    ontology.show()
    print("============================")
    name = 'xyl'
    result_list = ontology.fuzzy_search(search_string=name)
    for item in result_list:
        print(item)


#########################################################################
#### A very simple example of using this class
def peco_example():
    ontology = Ontology(filename='peco.obo',verbose=1)
    ontology.show()
    print("============================")
    name = 'light'
    result_list = ontology.fuzzy_search(search_string=name)
    for item in result_list:
        print(item)


#########################################################################
#### If class is run directly
def main():
    #psims_example()
    peco_example()

if __name__ == "__main__": main()
