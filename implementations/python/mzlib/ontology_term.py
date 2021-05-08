#!/usr/bin/env python3

#from __future__ import print_function
#import sys
#def eprint(*args, **kwargs):
#    print(*args, file=sys.stderr, **kwargs)

import logging
import re
import sys

#############################################################################
#### Ontology Term class
class OntologyTerm(object):


    #########################################################################
    #### Constructor
    def __init__(self, line_list=None, verbose=0):

        self.line_list = []
        self.verbose = verbose

        self.is_valid = False
        self.prefix = None
        self.identifier = None
        self.curie = None
        self.name = None
        self.value_type = None
        self.definition = None
        self.origin = None
        self.unparsable_line_list = []
        self.origin = None
        self.xref_list = []
        self.relationship_list = []
        self.parents = []
        self.children = []
        self.synonyms = []
        self.xrefs = []
        self.has_units = []
        self.has_order = None
        self.has_domain = None
        self.has_regexp = None
        self.is_obsolete = False
        self.namespaces = []
        self.subsets = []

        #### Mass modification-related properties
        self.monoisotopic_mass = None
        self.average_mass = None
        self.sites = {}
        self.extended_name = None

        self.n_errors = 0
        self.error_code = None
        self.error_message = None

        #### If we have been given an input line_list on construction, parse it right away
        if line_list is not None:
            self.parse(line_list=line_list)
        

    #########################################################################
    #### parse the line_list
    def parse(self, line_list=None, verbose=0):
        verboseprint = print if verbose>1 else lambda *a, **k: None
        if verbose > 1:
            logger = logging.getLogger()
            logger.setLevel(logging.DEBUG)

        #### Make sure there is a line_list to parse
        if line_list is not None:
            self.line_list = line_list

        #### Loop over the lines, processing each one
        for line in self.line_list:

            #### Just skip the [Term] line
            match = re.search("^\s*\[Term\]\s*$",line)
            if match: continue

            has_match = False

            #############################
            #### Process the id line
            match = re.search("^\s*id:",line)
            if match:
                match = re.search("^\s*id:\s*(\S+)\s*$",line)
                if match:
                    self.curie = match.groups()[0]
                    if ":" in self.curie:
                        self.prefix,self.identifier = self.curie.split(":",1)
                        has_match = True
                    else:
                        #logging.error("Unable to parse '%s'", line)
                        self.prefix = ''
                        self.identifier = self.curie
                        has_match = True
                else:
                    self.set_error("TermIdError",f"Unable to parse id line '{line}'")

            #############################
            #### Process the name line
            match = re.search("^\s*name:",line)
            if match:
                match = re.search("^\s*name:\s*(.+)\s*$",line)
                if match:
                    self.name = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermNameError",f"Unable to parse id line '{line}'")

            #############################
            #### Process the def line
            match = re.search("^\s*def:",line)
            if match:
                match = re.search("^\s*def:\s*(.+)\s*$",line)
                if match:
                    self.definition = match.groups()[0]
                    has_match = True
                    match = re.search('^"(.+)"\s*\[(.*)\].*$',self.definition)
                    if match:
                        self.definition = match.groups()[0]
                        self.origin = match.groups()[1]
                    else:
                        logging.error("Unable to parse definition string '%s'", self.definition)
                    has_match = True
                else:
                    self.set_error("TermDefinitionError",f"Unable to parse def line '{line}'")

            #############################
            #### Process the value-type line
            match = re.search("^\s*xref: value-type",line)
            if match:
                if self.value_type is not None:
                    logging.error("This term already has a type at line '%s'", line)
                match = re.search("^\s*xref: value-type:(\S+)",line)
                if match:
                    self.value_type = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermValueTypeError",f"Unable to parse value-type line '{line}'")

            #############################
            #### Process the binary-data-type line
            match = re.search("^\s*xref: binary-data-type",line)
            if match:
                match = re.search("^\s*xref: binary-data-type:(\S+)",line)
                if match:
                    pass
                    has_match = True
                else:
                    self.set_error("TermBinaryDataTypeError",f"Unable to parse binary-data-type line '{line}'")

            #############################
            #### Process the is_a line
            match = re.search("^\s*is_a",line)
            if match:
                match = re.search("^\s*is_a:\s*(\S+)",line)
                if match:
                    self.parents.append( { "type": "is_a", "curie": match.groups()[0] } )
                    has_match = True
                else:
                    self.set_error("TermIsAError",f"Unable to parse is_a line '{line}'")

            #############################
            #### Process the part_of line
            match = re.search("^relationship: part_of",line)
            if match:
                match = re.search("^\s*relationship:\s*part_of\s*(\S+)",line)
                if match:
                    self.parents.append( { "type": "part_of", "curie": match.groups()[0] } )
                    has_match = True
                else:
                    self.set_error("TermPartOfError",f"Unable to parse part_of line '{line}'")

            #############################
            #### Process the has_units line
            match = re.search("has_units",line)
            if match:
                match = re.search("^\s*relationship: has_units\s+(\S+)",line)
                if match:
                    self.has_units.append(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermHasUnitsError",f"Unable to parse has_units line '{line}'")

            #############################
            #### Process the has_order line
            match = re.search("has_order",line)
            if match:
                match = re.search("^\s*relationship: has_order\s+(\S+)",line)
                if match:
                    self.has_order = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermHasOrderError",f"Unable to parse has_order line '{line}'")

            #############################
            #### Process the has_domain line
            match = re.search("has_domain",line)
            if match:
                match = re.search("^\s*relationship: has_domain\s+(\S+)",line)
                if match:
                    self.has_domain = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermHasDomainError",f"Unable to parse has_domain line '{line}'")

            #############################
            #### Process the has_regexp line
            match = re.search("has_regexp",line)
            if match:
                match = re.search("^\s*relationship: has_regexp\s+(\S+)",line)
                if match:
                    self.has_regexp = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermHasRegexpError",f"Unable to parse has_regexp line '{line}'")

            #############################
            #### Process the is_obsolete line
            match = re.search("^\s*is_obsolete",line)
            if match:
                match = re.search("^\s*is_obsolete:\s*true\s*",line)
                if match:
                    self.is_obsolete = True
                    has_match = True
                else:
                    self.set_error("TermIsObsoleteError",f"Unable to parse is_obsolete line '{line}'")

            #############################
            #### Process the comment line
            match = re.search("^\s*comment",line)
            if match:
                match = re.search("^\s*comment:\s*(.+)\s*",line)
                if match:
                    self.comment = match.groups()[0]
                    has_match = True
                else:
                    self.set_error("TermCommentError",f"Unable to parse comment line '{line}'")

            #############################
            #### Process the synonym line
            match = re.search("^\s*synonym",line)
            if match:
                match = re.search('Japanese|Spanish',line)
                if match:
                    continue
                match = re.search('^\s*synonym:\s*"(.+)"\s*(\S+)\s*\[(.*)\].*$',line)
                if match:
                    self.synonyms.append( { "type": match.groups()[1], "term": match.groups()[0], "origin": match.groups()[2] } )
                    self.synonyms,match.groups()[0]
                    has_match = True
                else:
                    match = re.search('^\s*synonym:\s*"(.+)"\s*(\S+)\s+(\S+)\s*\[(.*)\]\s*$',line)
                    if match:
                        self.synonyms.append( { "type": match.groups()[1], "term": match.groups()[0], "domain": match.groups()[2], "origin": match.groups()[3] } )
                        self.synonyms,match.groups()[0]
                        has_match = True
                    else:
                        match = re.search('^\s*synonym:\s*"(.+)"\s*\[(.*)\]\s*$',line)
                        if match:
                            self.synonyms.append( { "type": 'unspecified', "term": match.groups()[0], "origin": match.groups()[1] } )
                            self.synonyms,match.groups()[0]
                            has_match = True
                        else:
                            #self.set_error("TermSynonymError",f"Unable to parse synonym line '{line}'")
                            logging.error("TermSynonymError, Unable to parse synonym line '%s'", line)
                            has_match = True

            #############################
            #### Process the alt_id line
            match = re.search("^\s*alt_id",line)
            if match:
                match = re.search('^\s*alt_id:\s*(\S+)\s*$',line)
                if match:
                    pass
                    has_match = True
                else:
                    self.set_error("TermAltIdError",f"Unable to parse alt_id line '{line}'")

            #############################
            #### Process the replaced_by line
            match = re.search("^\s*replaced_by",line)
            if match:
                match = re.search('^\s*replaced_by:\s*(\S+)\s*$',line)
                if match:
                    pass
                    has_match = True
                else:
                    self.set_error("TermReplacedByError",f"Unable to parse replaced_by line '{line}'")

            #############################
            #### Process the property_value line
            match = re.search("^\s*property_value",line)
            if match:
                match = re.search('^\s*property_value:\s*(.+)$',line)
                if match:
                    pass
                    has_match = True
                else:
                    self.set_error("TermPropertyValueError",f"Unable to parse property_value line '{line}'")

            #############################
            #### Process the other uninteresting lines
            match = re.search("^\s*(consider|disjoint_from|intersection_of|created_by|creation_date|equivalent_to|union_of)",line)
            if match:
                pass
                has_match = True

            #############################
            #### Process a namespace line
            match = re.search("^namespace:",line)
            if match:
                match = re.search("^\s*namespace:\s*(.+)\s*$",line)
                if match:
                    self.namespaces.append(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermNamespaceError",f"Unable to parse namespeace line '{line}'")

            #############################
            #### Process a subset line
            match = re.search("^subset:",line)
            if match:
                match = re.search("^\s*subset:\s*(.+)\s*$",line)
                if match:
                    self.subsets.append(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermSubsetError",f"Unable to parse subset line '{line}'")


            ####################################################################################
            #### Process mass modification data

            #### Parse xref: delta_mono_mass
            match = re.search(r"^xref:\s*delta_mono_mass",line)
            if has_match is False and match:
                match = re.search(r'^xref:\s*delta_mono_mass\s+\"\s*([\+\-\.\d]+)\s*\"',line)
                if match:
                    self.monoisotopic_mass = float(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermDeltaMonoMassError",f"Unable to parse xref line '{line}'")

            #### Parse xref: delta_avge_mass
            match = re.search(r"^xref:\s*delta_avge_mass",line)
            if has_match is False and match:
                match = re.search(r'^xref:\s*delta_avge_mass\s+\"\s*([\+\-\.\d]+)\s*\"',line)
                if match:
                    self.average_mass = float(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermDeltaAvgMassError",f"Unable to parse xref line '{line}'")

            #### Parse xref spec_NN_site
            match = re.search(r"^xref:\s+spec_\d+_site",line)
            if has_match is False and match:
                match = re.search(r'^xref:\s+spec_\d+_site\s+\"\s*(.+)\s*\"',line)
                if match:
                    self.sites[match.groups()[0]] = 1
                    has_match = True
                else:
                    self.set_error("TermSpecSiteError",f"Unable to parse xref line '{line}'")


            #############################
            #### Process an other kind of xref line
            match = re.search("^xref",line)
            if has_match is False and match:
                match = re.search("^\s*xref:\s+(\S+)\s*.*$",line)
                if match:
                    self.xrefs.append(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermXrefError",f"Unable to parse xref line '{line}'")

            #############################
            #### Process an other kind of relationship line
            match = re.search("^relationship",line)
            if has_match is False and match:
                match = re.search("^\s*relationship: (.+)\s*$",line)
                if match:
                    self.relationship_list.append(match.groups()[0])
                    has_match = True
                else:
                    self.set_error("TermRelationshipError",f"Unable to parse relationship line '{line}'")

            #############################
            #### If no match was found, add it to a pile of stuff we don't know how to deal with
            if has_match is False:
                self.unparsable_line_list.append(line)

        #### Set the is_valid state
        if self.n_errors == 0:
           self.is_valid = True

        else:
            self.is_valid = False
            logging.critical("Number of errors while parsing term '%s': %i", self.name, self.n_errors)
 
        if self.n_errors > 0 or len(self.unparsable_line_list) > 0:
            print("=====================")
            self.show()
            sys.exit()


    #########################################################################
    # Set the term to the error state
    def set_error(self,error_code,error_message):
        self.error_code = error_code
        self.error_message = error_message
        self.n_errors += 1
        logging.error("(%s): %s", error_code, error_message)


    #########################################################################
    # Print out some information about the term
    def show(self):
        print(f"curie: {self.curie}")
        print(f"is_valid: {self.is_valid}")
        print(f"prefix: {self.prefix}")
        print(f"identifier: {self.identifier}")
        print(f"name: {self.name}")
        print(f"definition: {self.definition}")
        print(f"origin: {self.origin}")
        print(f"value_type: {self.value_type}")
        print(f"parents: {self.parents}")
        print(f"children: {self.children}")
        print(f"synonyms: {self.synonyms}")
        print(f"has_units: {self.has_units}")
        print(f"is_obsolete: {self.is_obsolete}")

        if self.monoisotopic_mass is not None:
            print(f"=monoisotopic_mass: {self.monoisotopic_mass}")
            print(f"=average_mass: {self.average_mass}")
            print(f"=sites: {self.sites}")

        print(f"Number of unparsable lines: {len(self.unparsable_line_list)}")
        if len(self.unparsable_line_list) > 0:
            i_other_line = 0
            for line in self.unparsable_line_list:
                print("  >"+line)
                i_other_line += 1
                if i_other_line > 15: break
            if i_other_line > 15:
                print("  >...")




#########################################################################
#### A very simple example of using this class
def example():
    pass


#########################################################################
#### If class is run directly
def main():
    example()

if __name__ == "__main__": main()
