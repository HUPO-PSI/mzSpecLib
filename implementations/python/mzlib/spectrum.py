#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re
import json

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

        self.attributes = []
        self.attribute_dict = {}
        self.group_dict = {}

        self.peak_list = []
        self.group_counter = 1

        #### When ingesting a foreign format, temporarily store those foreign attributes here
        self.foreign_attributes = {}


    #### Get the next group identifier
    def get_next_group_identifier(self):
        next = self.group_counter
        self.group_counter += 1
        return(str(next))


    #### Add an attribute to the list and update the lookup tables
    def add_attribute(self, key, value, group_identifier=None):
        items = [ key, value ]
        if group_identifier is not None:
            items.append(group_identifier)
        self.attributes.append(items)
        index = len(self.attributes)-1

        #### If there is already one of these, add it to the lists in attribute_dict
        if key in self.attribute_dict:
            self.attribute_dict[key]["indexes"].append(index)
            if group_identifier is not None:
                self.attribute_dict[key]["groups"].append(group_identifier)

        #### Otherwise, create the entry in attribute_dict
        else:
            if group_identifier is not None:
                self.attribute_dict[key] = { "indexes": [ index ], "groups": [ group_identifier ] }
            else:
                self.attribute_dict[key] = { "indexes": [ index ], "groups": [ ] }

        #### If there is a group_identifier, then update the group_dict
        if group_identifier is not None:
            #### If this group already has one or more entries, add to it
            if group_identifier in self.group_dict:
                self.group_dict[group_identifier].append(index)
            #### Else create and entry for the group_identifier
            else:
                self.group_dict[group_identifier] = [ index ]

        return()


    #### Parse a list buffer of lines from a MSP-style spectrum entry, creating
    #### a dict of attributes and a list of peaks
    def parse(self, buffer, spectrum_index=None):

        #### Start in the header section of the entry
        in_header = True
        
        #### Reset all spectrum properties in case there is already data here
        self.foreign_attributes = {}
        self.peak_list = []

        #### Loop through each line in the buffered list
        for line in buffer:

            #print(line)

            #### If in the the header portion of the entry
            if in_header:

                #### Extract the key,value pair by splitting on the *first* colon with optional whitespace
                match = re.match("\s*#",line)
                if match:
                    continue
                elif line.count(":") > 0:
                    key, value = re.split(":\s*", line, 1)
                elif line.count("=") > 0:
                    key, value = re.split("=\s*", line, 1)
                elif line.count("\t") > 0:
                    print("ERROR: Looks like peaks in the header???")
                    in_header = False
                else:
                    key = line
                    value = None

                #### Adds the key-value pair to the dictionary
                self.foreign_attributes[key] = value
                
                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key == "Num peaks":
                    in_header = False
                    
                #### The "Comment" key requires special parsing
                if key == "Comment":

                    #### Remove it from attributes
                    del(self.foreign_attributes[key])

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
                            self.foreign_attributes[comment_key] = comment_value
                            #print(f"{comment_key}={comment_value}")
                        #### Otherwise just store the key with a null value
                        else:
                            self.foreign_attributes[item] = None
                            #print(f"{item}")

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                values = re.split(r'\s+',line)
                interpretations = ""
                if len(values) == 1:
                    mz = values
                    intensity = "1"
                if len(values) == 2:
                    mz, intensity = values
                elif len(values) == 3:
                    mz, intensity, interpretations = values
                elif len(values) > 3:
                    mz, intensity, interpretations = values[0:2]
                else:
                    mz = "1"
                    intensity = "1"

                interpretations = interpretations.strip('"')

                #### Add to the peak list
                self.peak_list.append( [ mz, intensity, interpretations ] )

        #### Now convert the foreign attributes to standard ones
        if spectrum_index is not None:
            self.add_attribute("MS:1008014|spectrum index", spectrum_index)
        self.convert_foreign_attributes()

        return(self)


    def convert_foreign_attributes(self):
        """
        convert_foreign_attributes - Convert the foreign attributes into standard ones

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Define the translation from keys to CV terms
        leader_terms = {
            "Name": "MS:1008013|spectrum name",
        }
        other_terms = {
            "MW": "MS:1008010|molecular mass",
            "ExactMass": "MS:1008010|molecular mass",
            "Charge": "MS:1000041|charge state",
            "Parent": "MS:1000744|selected ion m/z",
            "ObservedPrecursorMZ": "MS:1000744|selected ion m/z",
            "Single": [ "MS:1008015|spectrum aggregation type", "MS:1008016|singleton spectrum" ],
            "Consensus": [ "MS:1008015|spectrum aggregation type", "MS:1008017|consensus spectrum" ],
            "PrecursorMonoisoMZ": "MS:1008032|theoretical monoisotopic m/z",
            "Mz_exact": "MS:1008032|theoretical monoisotopic m/z",
            "Mz_av": "MS:1008033|theoretical average m/z",
            "Inst": { "it": [ [ "MS:1000044|dissociation method", "MS:1002472|trap-type collision-induced dissociation" ] ],
                      "hcd": [ [ "MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation" ] ] },
            "Pep": { "Tryptic": [ [ "MS:1008030|number of enzymatic termini", 2 ], [ "MS:1001045|cleavage agent name", "MS:1001251|Trypsin" ] ],
                     "N-Semitryptic": [ [ "MS:1008030|number of enzymatic termini", 1 ], [ "MS:1001045|cleavage agent name", "MS:1001251|Trypsin" ] ],
                     "C-Semitryptic": [ [ "MS:1008030|number of enzymatic termini", 1 ], [ "MS:1001045|cleavage agent name", "MS:1001251|Trypsin" ] ],
                     "Tryptic/miss_good_confirmed": [ [ "MS:1008030|number of enzymatic termini", 2 ],
                        [ "MS:1008034|number of missed cleavages", "0" ],
                        [ "MS:1001045|cleavage agent name", "MS:1001251|Trypsin" ] ],
                     "Tryptic/miss_bad_confirmed": [ [ "MS:1008030|number of enzymatic termini", 2 ],
                        [ "MS:1008034|number of missed cleavages", ">0" ],
                        [ "MS:1001045|cleavage agent name", "MS:1001251|Trypsin" ] ],
                    },
            "Spec": { "Consensus": [ [ "MS:1008015|spectrum aggregation type", "MS:1008017|consensus spectrum" ] ] },
            "Scan": "MS:1008035|scan number",
            "Origfile": "MS:1009008|source file",
            "Sample": "MS:1000002|sample name",
            "Filter": "MS:1000512|filter string",
            "FTResolution": "MS:1000028|detector resolution",
            "Protein": "MS:1000885|protein accession",
            "ms1PrecursorAb": "MS:1009010|previous MS1 scan precursor intensity",
            "Precursor1MaxAb": "MS:1009011|precursor apex intensity",
            "Purity": "MS:1009013|isolation window precursor purity",
            "Unassigned": "MS:1008027|top 20 peak unassigned intensity fraction",
            "Unassign_all": "MS:1008025|total unassigned intensity fraction",
            "Protein": "MS:1000885|protein accession",
            "Mods": "MS:1001471|peptide modification details",
            "BasePeak": "MS:1000505|base peak intensity",
            "Naa": "MS:1008027|number of residues",
            "Num peaks": "MS:1008040|number of peaks",
        }

        #### Add special terms that we want to start off with
        for term in leader_terms:
            if term in self.foreign_attributes:
                self.add_attribute(leader_terms[term], self.foreign_attributes[term])
            else:
                self.add_attribute("ERROR", f"Required term {leader_terms[term]} is missing")

        #### Translate the rest of the known attributes and collect unknown ones
        unknown_terms = []
        for attribute in self.foreign_attributes:

            #### Skip a leader term that we already processed
            if attribute in leader_terms: continue

            #### If this is in the list of generic terms, process it
            if attribute in other_terms:

                #### If the original attribute has no value, if mapping permits this, go ahead
                if self.foreign_attributes[attribute] is None:
                    if type(other_terms[attribute]) is list:
                        self.add_attribute(other_terms[attribute][0],other_terms[attribute][1])
                    else:
                        self.add_attribute("ERROR", f"Term {attribute} found without a value")
                        unknown_terms.append(attribute)

                #### If the term mapping value is an ordinary string, then just substitute the key
                elif type(other_terms[attribute]) is str:
                    self.add_attribute(other_terms[attribute], self.foreign_attributes[attribute])
                        
                #### Otherwise assume it is a dict of possible allowed values
                else:
                    #### If the value of the original attribute is in the list of allowed values
                    if self.foreign_attributes[attribute] in other_terms[attribute]:
                        #### If the mapping is a plain string, add it
                        if type(other_terms[attribute][self.foreign_attributes[attribute]]) is str:
                            self.add_attribute(other_terms[attribute][self.foreign_attributes[attribute]].split("="))
                        #### Or if it is a list, then there are multiple terms to add within a group
                        elif type(other_terms[attribute][self.foreign_attributes[attribute]]) is list:
                            if len(other_terms[attribute][self.foreign_attributes[attribute]]) == 1:
                                for item in other_terms[attribute][self.foreign_attributes[attribute]]:
                                    self.add_attribute(item[0],item[1])
                            else:
                                group_identifier = self.get_next_group_identifier()
                                for item in other_terms[attribute][self.foreign_attributes[attribute]]:
                                    self.add_attribute(item[0],item[1],group_identifier)
                        else:
                            self.add_attribute("ERROR", f"Internal error. Unexpected datatype in other_terms for {attribute}")
                            unknown_terms.append(attribute)
                    else:
                        self.add_attribute("ERROR", f"{self.foreign_attributes[attribute]} is not a defined option for {attribute}")
                        unknown_terms.append(attribute)

            #### Handle special logic cases

            #### Expand the HCD attribute
            elif attribute == "HCD":
                self.add_attribute("MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation")
                found_match = 0
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*ev", self.foreign_attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        found_match = 1
                        group_identifier = self.get_next_group_identifier()
                        self.add_attribute("MS:1000045|collision energy", match.group(1), group_identifier)
                        self.add_attribute("UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
                    match = re.match("([\d\.]+)\s*%", self.foreign_attributes[attribute])
                    if match is not None:
                        found_match = 1
                        group_identifier = self.get_next_group_identifier()
                        self.add_attribute("MS:1000045|collision energy", match.group(1), group_identifier)
                        self.add_attribute("UO:0000000|unit", "UO:0000187|percent", group_identifier)
                    if found_match == 0:
                        self.add_attribute("ERROR", f"Unable to parse self.foreign_attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Collision_energy attribute
            elif attribute == "Collision_energy":
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("([\d\.]+)", self.foreign_attributes[attribute])
                    if match is not None:
                        group_identifier = self.get_next_group_identifier()
                        self.add_attribute("MS:1000045|collision energy", match.group(1), group_identifier)
                        self.add_attribute("UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
                    else:
                        self.add_attribute("ERROR", f"Unable to parse self.foreign_attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the RT attribute
            elif attribute == "RT":
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*(\D*)", self.foreign_attributes[attribute])
                    if match is not None:
                        if match.group(2):
                            self.add_attribute("ERROR", f"Need more RT parsing code to handle this value")
                            unknown_terms.append(attribute)
                        else:
                            group_identifier = self.get_next_group_identifier()
                            self.add_attribute("MS:1000894|retention time", match.group(1), group_identifier)
                            #### If the value is greater than 250, assume it must be seconds
                            if float(match.group(1)) > 250:
                                self.add_attribute("UO:0000000|unit", "UO:0000010|second", group_identifier)
                            #### Although normally assume minutes
                            else:
                                self.add_attribute("UO:0000000|unit", "UO:0000031|minute", group_identifier)
                    else:
                        self.add_attribute("ERROR", f"Unable to parse self.foreign_attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the ms2IsolationWidth attribute
            elif attribute == "ms2IsolationWidth":
                if self.foreign_attributes[attribute] is not None:
                    group_identifier = self.get_next_group_identifier()
                    self.add_attribute("MS:1000828|isolation window lower offset", str(float(self.foreign_attributes[attribute])/2), group_identifier)
                    self.add_attribute("UO:0000000|unit", "MS:1000040|m/z", group_identifier)
                    group_identifier = self.get_next_group_identifier()
                    self.add_attribute("MS:1000829|isolation window upper offset", str(float(self.foreign_attributes[attribute])/2), group_identifier)
                    self.add_attribute("UO:0000000|unit", "MS:1000040|m/z", group_identifier)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Mz_diff attribute
            elif attribute == "Mz_diff":
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("([\-\+e\d\.]+)\s*ppm", self.foreign_attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        group_identifier = self.get_next_group_identifier()
                        self.add_attribute("MS:1001975|delta m/z", match.group(1), group_identifier)
                        self.add_attribute("UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
                    else:
                        match = re.match("([\-\+e\d\.]+)\s*", self.foreign_attributes[attribute])
                        if match is not None:
                            group_identifier = self.get_next_group_identifier()
                            self.add_attribute("MS:1001975|delta m/z", match.group(1), group_identifier)
                            self.add_attribute("UO:0000000|unit", "MS:1000040|m/z", group_identifier)
                        else:
                            self.add_attribute("ERROR", f"Unable to parse {self.foreign_attributes[attribute]} in {attribute}")
                            unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Dev_ppm attribute
            elif attribute == "Dev_ppm":
                if self.foreign_attributes[attribute] is not None:
                    group_identifier = self.get_next_group_identifier()
                    self.add_attribute("MS:1001975|delta m/z", self.foreign_attributes[attribute], group_identifier)
                    self.add_attribute("UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Fullname attribute
            elif attribute == "Fullname":
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("([A-Z\-\*])\.([A-Z]+)\.([A-Z\-\*])/*([\d]*)", self.foreign_attributes[attribute])
                    if match is not None:
                        self.add_attribute("MS:1000888|unmodified peptide sequence", match.group(2))
                        self.add_attribute("MS:1001112|n-terminal flanking residue", match.group(1))
                        self.add_attribute("MS:1001113|c-terminal flanking residue", match.group(3))
                        if match.group(4):
                            self.add_attribute("MS:1000041|charge state", match.group(4))
                    else:
                        self.add_attribute("ERROR", f"Unable to parse {self.foreign_attributes[attribute]} in {attribute} at E2355")
                        unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Nrep attribute
            elif attribute == "Nrep" or attribute == "Nreps":
                if self.foreign_attributes[attribute] is not None:
                    match = re.match("(\d+)/(\d+)", self.foreign_attributes[attribute])
                    if match is not None:
                        self.add_attribute("MS:1009020|number of replicate spectra used", match.group(1))
                        self.add_attribute("MS:1009021|number of replicate spectra available", match.group(2))
                    else:
                        match = re.match("(\d+)", self.foreign_attributes[attribute])
                        if match is not None:
                            self.add_attribute("MS:1008020|number of replicate spectra used", match.group(1))
                            self.add_attribute("MS:1008019|number of replicate spectra available", match.group(1))
                        else:
                            self.add_attribute("ERROR", f"Unable to parse {self.foreign_attributes[attribute]} in {attribute} at E2455")
                            unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Fullname attribute
            elif attribute == "Organism":
                if self.foreign_attributes[attribute] is not None:
                    value = self.foreign_attributes[attribute]
                    value = value.strip('"')
                    species_map = {
                        "human": [ [ "MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9606|Homo sapiens" ],
                            ["MS:1001469|taxonomy: scientific name", "Homo sapiens" ],
                            ["MS:1001468|taxonomy: common name", "human" ] ],
                        "zebrafish": [ [ "MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:7955|Danio rerio" ],
                            ["MS:1001469|taxonomy: scientific name", "Danio rerio" ],
                            ["MS:1001468|taxonomy: common name", "zebra fish" ] ],
                        "chicken": [ [ "MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9031|Gallus gallus" ],
                            ["MS:1001469|taxonomy: scientific name", "Gallus gallus" ],
                            ["MS:1001468|taxonomy: common name", "chicken" ] ],
                        }
                    if value in species_map:
                        group_identifier = self.get_next_group_identifier()
                        for item in species_map[value]:
                            self.add_attribute(item[0],item[1],group_identifier)
      
                    else:
                        self.add_attribute("ERROR", f"Unable to parse {self.foreign_attributes[attribute]} in {attribute} at E2355")
                        unknown_terms.append(attribute)
                else:
                    self.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Otherwise add this term to the list of attributes that we don't know how to handle
            else:
                unknown_terms.append(attribute)


        #### Perform some cleanup for poorly annotated libraries
        if "MS:1000888|unmodified peptide sequence" not in self.attribute_dict:
            if "MS:1008013|spectrum name" in self.attribute_dict:
                lookup = self.attribute_dict["MS:1008013|spectrum name"]
                name = self.attributes[ lookup["indexes"][0] ][1]
                match = re.match("(.+)/(\d+)",name)
                if match:
                    self.add_attribute("MS:1000888|unmodified peptide sequence", match.group(1))
                    self.add_attribute("MS:1000041|charge state", match.group(2))

        #### Handle the uninterpretable terms
        for attribute in unknown_terms:
            if self.foreign_attributes[attribute] is None:
                self.add_attribute("MS:1009900|other attribute name", attribute)
            else:
                group_identifier = self.get_next_group_identifier()
                self.add_attribute("MS:1009900|other attribute name", attribute, group_identifier)
                self.add_attribute("MS:1009902|other attribute value", self.foreign_attributes[attribute], group_identifier)

        return()


    def write(self, format="text"):
        """
        write - Write out the spectrum in any of the supported formats

        Extended description of function.

        Parameters
        ----------

        Returns
        -------
        int
            Description of return value
        """

        #### Set a buffer to fill with string data
        buffer = ''

        #### Make the format string lower case to facilitate comparisons
        format = format.lower()

        #### If the format is text
        if format == "text":
            for attribute in self.attributes:
                if len(attribute) == 2:
                    buffer += f"{attribute[0]}={attribute[1]}\n"
                elif len(attribute) == 3:
                    buffer += f"[{attribute[2]}]{attribute[0]}={attribute[1]}\n"
                else:
                    print("-ERROR: attribute has wrong number of elements:")
                    print(attribute)
            for peak in self.peak_list:
                buffer += "\t".join(peak)+"\n"

        #### If the format is TSV
        elif format == "tsv" or format == "csv":

            #### Set the appropriate delimiter for format
            delimiter = "\t"
            if format == "csv": delimiter = ","

            #### Create the header line and columns line
            buffer += "# Spectrum attributes\n"
            buffer += "cv_param_group\taccession\tname\tvalue_accession\tvalue\n"

            for attribute in self.attributes:
                if len(attribute) == 2:
                    key,value = attribute
                    cv_param_group = ''
                elif len(attribute) == 3:
                    key,value,cv_param_group = attribute
                    if format == "csv" and ',' in str(value):
                        value = '"' + value + '"'
                else:
                    print("ERROR: Unsupported number of items in attribute")
                    print(attribute)
                    sys.exit(10)
                components = key.split('|',1)
                if len(components) == 2:
                    accession,name = components
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    sys.exit(10)
                components = str(value).split('|',1)
                if len(components) == 2:
                    value_accession,value = components
                    value = str(value)
                    if format == "csv" and ',' in value:
                        value = '"' + value + '"'
                elif len(components) == 1:
                    value = str(value)
                    if format == "csv" and ',' in value:
                        value = '"' + value + '"'
                    value_accession = ''
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    sys.exit(10)

                #### Create the data line
                buffer += delimiter.join([cv_param_group,accession,name,value_accession,value])+"\n"

            #### Create the header line and columns line
            buffer += "# Peak list\n"
            buffer += "mz\tintensity\tinterpretation\n"

            #### Write out the peak list
            for peak in self.peak_list:
                mz,intensity,interpretation = peak
                if format == "csv" and ',' in str(interpretation):
                    interpretation = '"' + interpretation + '"'
                buffer += delimiter.join([mz,intensity,interpretation])+"\n"

        #### If the format is JSON
        elif format == "json":
            mzs = []
            intensities = []
            interpretations = []
            for peak in self.peak_list:
                mzs.append(float(peak[0]))
                intensities.append(float(peak[1]))
                interpretations.append(peak[2])

            #### Organize the attributes from the simple list into the appropriate JSON format
            attributes = []
            for attribute in self.attributes:
                reformed_attribute = {}
                if len(attribute) == 2:
                    key,value = attribute
                elif len(attribute) == 3:
                    key,value,cv_param_group = attribute
                    reformed_attribute['cv_param_group'] = cv_param_group
                else:
                    print("ERROR: Unsupported number of items in attribute")
                    print(attribute)
                    sys.exit(10)
                components = key.split('|',1)
                if len(components) == 2:
                    accession,name = components
                    reformed_attribute['accession'] = accession
                    reformed_attribute['name'] = name
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    sys.exit(10)
                components = str(value).split('|',1)
                if len(components) == 2:
                    value_accession,value = components
                    reformed_attribute['value_accession'] = value_accession
                    reformed_attribute['value'] = value
                elif len(components) == 1:
                    reformed_attribute['value'] = value
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    sys.exit(10)
                attributes.append(reformed_attribute)

            spectrum = { "attributes": attributes, "mzs": mzs, "intensities": intensities,
                "interpretations": interpretations }
            buffer = json.dumps(spectrum,sort_keys=True,indent=2)

        #### Otherwise we don't know this format
        else:
            raise ValueError(f"ERROR: Unrecogized format '{format}'")

        return(buffer)






