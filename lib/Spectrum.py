#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re

debug = True

#A class that holds data for each spectrum that is read from the SpectralLibrary class
class Spectrum:

    #### Constructor
    def __init__(self):
        """
        __init__ - SpectrumLibrary constructor

        Parameters
        ----------

        """

        self.attributes = {}
        self.peak_list = []
        self.group_counter = 1


    #### Get the next group identifier
    def get_next_group_identifier(self):
        next = self.group_counter
        self.group_counter += 1
        return(str(next))

    #### Parse a list buffer of lines from a MSP-style spectrum entry, creating
    #### a dict of attributes and a list of peaks
    def parse(self, buffer):

        #### Start in the header section of the entry
        in_header = True
        
        #### Reset all spectrum properties in case there is already data here
        self.attributes = {}
        self.peak_list = []

        #### Loop through each line in the buffered list
        for line in buffer:

            #### If in the the header portion of the entry
            if in_header:

                #### Extract the key,value pair by splitting on the *first* colon with optional whitespace
                key, value = re.split(":\s*", line, 1)
                
                #### Adds the key-value pair to the dictionary
                self.attributes[key] = value
                
                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key == "Num peaks":
                    in_header = False
                    
                #### The "Comment" key requires special parsing
                if key == "Comment":

                    #### Remove it from attributes
                    del(self.attributes[key])

                    #### Split on all spaces
                    comment_items = re.split(" ", value)

                    #### Any spaces within quotes are then de-split
                    fixed_comment_items = []
                    quote_counter = 0
                    new_item = ""
                    for item in comment_items:
                        if new_item > "":
                            new_item = new_item + " "
                        new_item = new_item + item
                        n_quotes = new_item.count('"')
                        if n_quotes/2 == int(n_quotes/2):
                            fixed_comment_items.append(new_item)
                            new_item = ""

                    #### Try to split each item on the first = character and store in attributes
                    for item in fixed_comment_items:
                        #### If there is an = then split on the first one and store key and value
                        if item.count("=") > 0:
                            comment_key, comment_value = item.split("=", 1)
                            self.attributes[comment_key] = comment_value
                        #### Otherwise just store the key with a null value
                        else:
                            self.attributes[item] = None

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                mz, intensity, interpretations = line.split( "\t" )
                
                #### Add to the peak list
                self.peak_list.append( [ mz, intensity, interpretations ] )

        return(self)


    def write_as_SpecML(self):
        """
        write - Write the spectrum in long form SpecML text format

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Create an empty buffer
        buffer = []
        #### Define the translation from keys to CV terms
        leader_terms = {
            "Name": "MS:1009004|spectrum label",
            "MW": "MS:1009004|molecular weight",
            "Charge": "MS:1000041|charge state"
        }
        other_terms = {
            "Parent": "MS:1000744|selected ion m/z",
            "Single": "MS:1009030|representative spectrum type=MS:1009031|individual spectrum",
            "Consensus": "MS:1009030|representative spectrum type=MS:1009032|consensus spectrum",
            "PrecursorMonoisoMZ": "MS:1009009|theoretical m/z",
            "Inst": { "it": "MS:1000044|dissociation method=MS:1002472|trap-type collision-induced dissociation" },
            "Pep": { "Tryptic": [ "MS:1009040|number of enzymatic termini=2", "MS:1001045|cleavage agent name=MS:1001251|Trypsin" ] },
            "Scan": "MS:1009007|scan number",
            "Origfile": "MS:1009008|source file",
            "Sample": "MS:1009009|sample name",
            "Filter": "MS:1000512|filter string",
            "FTResolution": "MS:1000028|detector resolution",
            "Protein": "MS:1000885|protein accession",
            "ms1PrecursorAb": "MS:1009010|previous scan precursor intensity",
            "Precursor1MaxAb": "MS:1009011|precursor apex intensity",
            "Purity": "MS:1009013|isolation window precursor purity",
            "Unassigned": "MS:1009014|top 20 peak unassigned intensity fraction",
            "Unassign_all": "MS:1009015|total unassigned intensity fraction",
            "Protein": "MS:1000885|protein accession",
            "Mods": "MS:1001471|peptide modification details",
            "Num peaks": "MS:1009006|number of peaks",
        }

        #### Start appending data
        buffer.append("MS:1009003|spectrum index=1")

        #### Add two special fields that we want to start off with
        for term in leader_terms:
            if self.attributes[term]:
                buffer.append(leader_terms[term] + "=" + self.attributes[term])
            else:
                buffer.append(leader_terms[term] + "=" + "Missing!")

        #### Translate the rest of the known attributes and collect unknown ones
        unknown_terms = []
        for attribute in self.attributes:
            if attribute in leader_terms: continue
            if attribute in other_terms:

                #### If the original attribute has no value, if mapping permits this, go ahead
                if self.attributes[attribute] is None:
                    if other_terms[attribute].count("=") > 0:
                        buffer.append(other_terms[attribute])
                    else:
                        buffer.append("? " + attribute + ": ERROR value not in lookup")

                #### If the term mapping is an ordinary string
                elif type(other_terms[attribute]) is str:
                    buffer.append(other_terms[attribute] + "=" + self.attributes[attribute])
                        
                else:
                    if self.attributes[attribute] in other_terms[attribute]:
                        if type(other_terms[attribute][self.attributes[attribute]]) is str:
                            buffer.append(other_terms[attribute][self.attributes[attribute]])
                        elif type(other_terms[attribute][self.attributes[attribute]]) is list:
                            group_indentifier = self.get_next_group_identifier()
                            gid = f"[{group_indentifier}]"
                            for item in other_terms[attribute][self.attributes[attribute]]:
                                buffer.append(f"{gid}{item}")
                        else:
                            buffer.append("? ERROR: Unable to process attribute with E455")
                    else:
                        buffer.append("? " + attribute + "=" + self.attributes[attribute] + ": ERROR value not in lookup")
            #### Handle special logic cases
            elif attribute == "HCD":
                buffer.append("MS:1000044|dissociation method=MS:1000422|beam-type collision-induced dissociation")
                if self.attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*ev", self.attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        group_indentifier = self.get_next_group_identifier()
                        gid = f"[{group_indentifier}]"
                        buffer.append(f"{gid}MS:1000045|collision energy" + "=" + match.group(1))
                        buffer.append(f"{gid}UO:0000000|unit=UO:0000266|electronvolt")
            elif attribute == "ms2IsolationWidth":
                if self.attributes[attribute] is not None:
                    group_indentifier = self.get_next_group_identifier()
                    gid = f"[{group_indentifier}]"
                    buffer.append(f"{gid}MS:1000828|isolation window lower offset=" + str(float(self.attributes[attribute])/2))
                    buffer.append(f"{gid}MS:1000040|m/z")
                    group_indentifier = self.get_next_group_identifier()
                    gid = f"[{group_indentifier}]"
                    buffer.append(f"{gid}MS:1000829|isolation window upper offset=" + str(float(self.attributes[attribute])/2))
                    buffer.append(f"{gid}MS:1000040|m/z")

            elif attribute == "Mz_diff":
                if self.attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*ppm", self.attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        group_indentifier = self.get_next_group_identifier()
                        gid = f"[{group_indentifier}]"
                        buffer.append(f"{gid}MS:1001975|delta m/z" + "=" + match.group(1))
                        buffer.append(f"{gid}UO:0000000|unit=UO:0000169|parts per million")
                    else:
                        buffer.append("? ERROR: Cannot parse" + attribute + self.attributes[attribute])
                else:
                    buffer.append("? ERROR: Cannot parse" + attribute)

            elif attribute == "Fullname":
                if self.attributes[attribute] is not None:
                    match = re.match("([A-Z\.\*])\.([A-Z]+)\.([A-Z\.\*])", self.attributes[attribute])
                    if match is not None:
                        buffer.append("MS:1000888|unmodified peptide sequence" + "=" + match.group(2))
                        buffer.append("MS:1001112|n-terminal flanking residue" + "=" + match.group(1))
                        buffer.append("MS:1001113|c-terminal flanking residue" + "=" + match.group(3))
                    else:
                        buffer.append("? ERROR: Cannot parse" + attribute + self.attributes[attribute])
                else:
                    buffer.append("? ERROR: Cannot parse" + attribute)

            elif attribute == "Nrep":
                if self.attributes[attribute] is not None:
                    match = re.match("(\d+)/(\d+)", self.attributes[attribute])
                    if match is not None:
                        buffer.append("MS:1009020|number of replicate spectra used" + "=" + match.group(1))
                        buffer.append("MS:1009021|number of replicate spectra available" + "=" + match.group(2))
                    else:
                        buffer.append("? ERROR: Cannot parse" + attribute + self.attributes[attribute])
                else:
                    buffer.append("? ERROR: Cannot parse" + attribute)

            else:
                unknown_terms.append(attribute)

        for attribute in unknown_terms:
            if self.attributes[attribute] is None:
                buffer.append("OtherAttribute=" + attribute)
            else:
                buffer.append("OtherAttribute=" + attribute + "=" + self.attributes[attribute])

        return(buffer)







