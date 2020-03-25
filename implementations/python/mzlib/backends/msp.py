import re

from .base import SpectralLibraryBackendBase

leader_terms = {
    "Name": "MS:1003061|spectrum name",
}

other_terms = {
    "MW": "MS:1000224|molecular mass",
    "ExactMass": "MS:1000224|molecular mass",
    "Charge": "MS:1000041|charge state",
    "Parent": "MS:1000744|selected ion m/z",
    "ObservedPrecursorMZ": "MS:1000744|selected ion m/z",
    "Single": ["MS:1003065|spectrum aggregation type", "MS:1003066|singleton spectrum"],
    "Consensus": ["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"],
    "PrecursorMonoisoMZ": "MS:1003053|theoretical monoisotopic m/z",
    "Mz_exact": "MS:1003053|theoretical monoisotopic m/z",
    "Mz_av": "MS:1003054|theoretical average m/z",
    "Inst": {"it": [["MS:1000044|dissociation method", "MS:1002472|trap-type collision-induced dissociation"]],
             "hcd": [["MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation"]]},
    "Pep": {"Tryptic": [["MS:1003048|number of enzymatic termini", 2], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
            "N-Semitryptic": [["MS:1003048|number of enzymatic termini", 1], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
            "C-Semitryptic": [["MS:1003048|number of enzymatic termini", 1], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
            "Tryptic/miss_good_confirmed": [["MS:1003048|number of enzymatic termini", 2],
                                            ["MS:1003044|number of missed cleavages", "0"],
                                            ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
            "Tryptic/miss_bad_confirmed": [["MS:1003048|number of enzymatic termini", 2],
                                           ["MS:1003044|number of missed cleavages", ">0"],
                                           ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
            },
    "Spec": {"Consensus": [["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"]]},
    "Scan": "MS:1003057|scan number",
            "Origfile": "MS:1009008|source file",
            "Sample": "MS:1000002|sample name",
            "Filter": "MS:1000512|filter string",
            "FTResolution": "MS:1000028|detector resolution",
            "Protein": "MS:1000885|protein accession",
            "ms1PrecursorAb": "MS:1009010|previous MS1 scan precursor intensity",
            "Precursor1MaxAb": "MS:1009011|precursor apex intensity",
            "Purity": "MS:1009013|isolation window precursor purity",
            "Unassigned": "MS:1003080|top 20 peak unassigned intensity fraction",
            "Unassign_all": "MS:1003079|total unassigned intensity fraction",
            "Mods": "MS:1001471|peptide modification details",
            "BasePeak": "MS:1000505|base peak intensity",
            "Naa": "MS:1003043|number of residues",
            "Num peaks": "MS:1003059|number of peaks",
}

species_map = {
    "human": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9606|Homo sapiens"],
                ["MS:1001469|taxonomy: scientific name",
                    "Homo sapiens"],
                ["MS:1001468|taxonomy: common name", "human"]],
    "zebrafish": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:7955|Danio rerio"],
                    ["MS:1001469|taxonomy: scientific name",
                        "Danio rerio"],
                    ["MS:1001468|taxonomy: common name", "zebra fish"]],
    "chicken": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9031|Gallus gallus"],
                ["MS:1001469|taxonomy: scientific name",
                    "Gallus gallus"],
                ["MS:1001468|taxonomy: common name", "chicken"]],
}

class MSPSpectralLibrary(SpectralLibraryBackendBase):

    def __init__(self, filename):
        super(MSPSpectralLibrary, self).__init__(filename)
        self.index = None

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
        # if debug:
        #     eprint(f'INFO: Reading spectra from {filename}')
        # if filename is None:
        #     eprint("ERROR: Unable to read library with no filename")
        #     return(False)

        #### If an index hasn't been opened, open it now
        # if create_index is not None:
        #     if self.index is None:
        #         self.index = SpectrumLibraryIndex(
        #             library_filename=self.filename)
        #         self.index.create_index()

        #### Determine the filesize
        file_size = os.path.getsize(filename)
        # if debug:
        #     eprint(f"INFO: File size is {file_size}")

        with open(filename, 'r') as infile:
            state = 'header'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''

            # Required for counting file_offset manually (LF vs CRLF)
            infile.readline()
            file_offset_line_ending = len(infile.newlines) - 1
            infile.seek(0)

            # if debug:
            #     eprint("INFO: Reading..", end='', flush=True)
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                # file_offset = infile.tell()
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':
                    if re.match('Name: ', line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('Name: ', line):
                        if len(spectrum_buffer) > 0:
                            if create_index is not None:
                                self.index.add_spectrum(
                                    number=n_spectra + start_index, offset=spectrum_file_offset, name=spectrum_name, peptide_sequence=None)
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            # if int(n_spectra/1000) == n_spectra/1000:
                            #     self.index.commit()
                            #     percent_done = int(
                            #         file_offset/file_size*100+0.5)
                            #     eprint(str(percent_done)+"%..",
                            #            end='', flush=True)

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match('Name:\s+(.+)', line).group(1)

                    spectrum_buffer.append(line)

            # if create_index is not None:
            #     self.index.add_spectrum(
            #         number=n_spectra + start_index, offset=spectrum_file_offset, name=spectrum_name, peptide_sequence=None)
            #     self.index.commit()
            n_spectra += 1

            #### Flush the index
            # if create_index:
            #     self.index.commit()

        return(n_spectra)

    def _get_lines_for(self, offset):
        with open(self.filename, 'r') as infile:
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
                    if re.match('Name: ', line):
                        if len(spectrum_buffer) > 0:
                            return spectrum_buffer
                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match('Name:\s+(.+)', line).group(1)
                    spectrum_buffer.append(line)

            #### We will end up here if this is the last spectrum in the file
            return spectrum_buffer

    def parse(self, buffer, spectrum_index=None):

        #### Start in the header section of the entry
        in_header = True

        #### Reset all spectrum properties in case there is already data here
        attributes = {}
        peak_list = []

        #### Loop through each line in the buffered list
        for line in buffer:

            #print(line)

            #### If in the the header portion of the entry
            if in_header:

                #### Extract the key,value pair by splitting on the *first* colon with optional whitespace
                match = re.match("\s*#", line)
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
                attributes[key] = value

                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key == "Num peaks":
                    in_header = False

                #### The "Comment" key requires special parsing
                if key == "Comment":

                    #### Remove it from attributes
                    del attributes[key]
                    self._parse_comment(value, attributes)

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                values = re.split(r'\s+', line)
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
                peak_list.append([mz, intensity, interpretations])

        #### Now convert the format attributes to standard ones
        if spectrum_index is not None:
            self.add_attribute("MS:1003062|spectrum index", spectrum_index)
        spectrum = self.make_spectrum(peak_list, attributes)

        return spectrum

    def _parse_comment(self, value, attributes):
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
                attributes[comment_key] = comment_value
                #print(f"{comment_key}={comment_value}")
            #### Otherwise just store the key with a null value
            else:
                attributes[item] = None

    def make_spectrum(self, peak_list, attributes):
        spectrum = self._new_spectrum()
        spectrum.peak_list = peak_list

        group_counter = self._make_counter()

        #### Add special terms that we want to start off with
        for term in leader_terms:
            if term in attributes:
                spectrum.add_attribute(
                    leader_terms[term], attributes[term])
            else:
                spectrum.add_attribute(
                    "ERROR", f"Required term {leader_terms[term]} is missing")
        #### Translate the rest of the known attributes and collect unknown ones
        unknown_terms = []
        for attribute in attributes:

            #### Skip a leader term that we already processed
            if attribute in leader_terms:
                continue

            #### If this is in the list of generic terms, process it
            if attribute in other_terms:
                                #### If the original attribute has no value, if mapping permits this, go ahead
                if attributes[attribute] is None:
                    if type(other_terms[attribute]) is list:
                        spectrum.add_attribute(
                            other_terms[attribute][0], other_terms[attribute][1])
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Term {attribute} found without a value")
                        unknown_terms.append(attribute)

                #### If the term mapping value is an ordinary string, then just substitute the key
                elif type(other_terms[attribute]) is str:
                    spectrum.add_attribute(
                        other_terms[attribute], attributes[attribute])

                #### Otherwise assume it is a dict of possible allowed values
                else:
                    #### If the value of the original attribute is in the list of allowed values
                    if attributes[attribute] in other_terms[attribute]:
                        #### If the mapping is a plain string, add it
                        if type(other_terms[attribute][attributes[attribute]]) is str:
                            spectrum.add_attribute(
                                other_terms[attribute][attributes[attribute]].split("="))
                        #### Or if it is a list, then there are multiple terms to add within a group
                        elif type(other_terms[attribute][attributes[attribute]]) is list:
                            if len(other_terms[attribute][attributes[attribute]]) == 1:
                                for item in other_terms[attribute][attributes[attribute]]:
                                    spectrum.add_attribute(item[0], item[1])
                            else:
                                group_identifier = spectrum.get_next_group_identifier()
                                for item in other_terms[attribute][attributes[attribute]]:
                                    spectrum.add_attribute(
                                        item[0], item[1], group_identifier)
                        else:
                            spectrum.add_attribute(
                                "ERROR", f"Internal error. Unexpected datatype in other_terms for {attribute}")
                            unknown_terms.append(attribute)
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"{attributes[attribute]} is not a defined option for {attribute}")
                        unknown_terms.append(attribute)
            elif attribute == "HCD":
                spectrum.add_attribute("MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation")
                found_match = 0
                if attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*ev", attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        found_match = 1
                        group_identifier = next(group_counter)
                        spectrum.add_attribute("MS:1000045|collision energy", match.group(1), group_identifier)
                        spectrum.add_attribute("UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
                    match = re.match("([\d\.]+)\s*%", attributes[attribute])
                    if match is not None:
                        found_match = 1
                        group_identifier = next(group_counter)
                        spectrum.add_attribute("MS:1000045|collision energy", match.group(1), group_identifier)
                        spectrum.add_attribute("UO:0000000|unit", "UO:0000187|percent", group_identifier)
                    if found_match == 0:
                        spectrum.add_attribute("ERROR", f"Unable to parse attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute("ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)
            elif attribute == "Collision_energy":
                if attributes[attribute] is not None:
                    match = re.match(
                        "([\d\.]+)", attributes[attribute])
                    if match is not None:
                        group_identifier = next(group_counter)
                        spectrum.add_attribute(
                            "MS:1000045|collision energy", match.group(1), group_identifier)
                        spectrum.add_attribute(
                            "UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Unable to parse attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)
            elif attribute == "RT":
                if attributes[attribute] is not None:
                    match = re.match("([\d\.]+)\s*(\D*)",
                                     attributes[attribute])
                    if match is not None:
                        if match.group(2):
                            spectrum.add_attribute(
                                "ERROR", f"Need more RT parsing code to handle this value")
                            unknown_terms.append(attribute)
                        else:
                            group_identifier = next(group_counter)
                            spectrum.add_attribute(
                                "MS:1000894|retention time", match.group(1), group_identifier)
                            #### If the value is greater than 250, assume it must be seconds
                            if float(match.group(1)) > 250:
                                spectrum.add_attribute(
                                    "UO:0000000|unit", "UO:0000010|second", group_identifier)
                            #### Although normally assume minutes
                            else:
                                spectrum.add_attribute(
                                    "UO:0000000|unit", "UO:0000031|minute", group_identifier)
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Unable to parse attributes[attribute] in attribute")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            elif attribute == "ms2IsolationWidth":
                if attributes[attribute] is not None:
                    group_identifier = next(group_counter)
                    spectrum.add_attribute("MS:1000828|isolation window lower offset", str(
                        float(attributes[attribute])/2), group_identifier)
                    spectrum.add_attribute("UO:0000000|unit",
                                       "MS:1000040|m/z", group_identifier)
                    group_identifier = next(group_counter)
                    spectrum.add_attribute("MS:1000829|isolation window upper offset", str(
                        float(attributes[attribute])/2), group_identifier)
                    spectrum.add_attribute("UO:0000000|unit",
                                       "MS:1000040|m/z", group_identifier)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Mz_diff attribute
            elif attribute == "Mz_diff":
                if attributes[attribute] is not None:
                    match = re.match(
                        "([\-\+e\d\.]+)\s*ppm", attributes[attribute], flags=re.IGNORECASE)
                    if match is not None:
                        group_identifier = next(group_counter)
                        spectrum.add_attribute(
                            "MS:1001975|delta m/z", match.group(1), group_identifier)
                        spectrum.add_attribute(
                            "UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
                    else:
                        match = re.match(
                            "([\-\+e\d\.]+)\s*", attributes[attribute])
                        if match is not None:
                            group_identifier = next(group_counter)
                            spectrum.add_attribute(
                                "MS:1001975|delta m/z", match.group(1), group_identifier)
                            spectrum.add_attribute(
                                "UO:0000000|unit", "MS:1000040|m/z", group_identifier)
                        else:
                            spectrum.add_attribute(
                                "ERROR", f"Unable to parse {attributes[attribute]} in {attribute}")
                            unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Dev_ppm attribute
            elif attribute == "Dev_ppm":
                if attributes[attribute] is not None:
                    group_identifier = next(group_counter)
                    spectrum.add_attribute(
                        "MS:1001975|delta m/z", attributes[attribute], group_identifier)
                    spectrum.add_attribute(
                        "UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Fullname attribute
            elif attribute == "Fullname":
                if attributes[attribute] is not None:
                    match = re.match(
                        "([A-Z\-\*])\.([A-Z]+)\.([A-Z\-\*])/*([\d]*)", attributes[attribute])
                    if match is not None:
                        spectrum.add_attribute(
                            "MS:1000888|unmodified peptide sequence", match.group(2))
                        spectrum.add_attribute(
                            "MS:1001112|n-terminal flanking residue", match.group(1))
                        spectrum.add_attribute(
                            "MS:1001113|c-terminal flanking residue", match.group(3))
                        if match.group(4):
                            spectrum.add_attribute(
                                "MS:1000041|charge state", match.group(4))
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Unable to parse {attributes[attribute]} in {attribute} at E2355")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Expand the Nrep attribute
            elif attribute == "Nrep" or attribute == "Nreps":
                if attributes[attribute] is not None:
                    match = re.match(
                        "(\d+)/(\d+)", attributes[attribute])
                    if match is not None:
                        spectrum.add_attribute(
                            "MS:1009020|number of replicate spectra used", match.group(1))
                        spectrum.add_attribute(
                            "MS:1009021|number of replicate spectra available", match.group(2))
                    else:
                        match = re.match(
                            "(\d+)", attributes[attribute])
                        if match is not None:
                            spectrum.add_attribute(
                                "MS:1003070|number of replicate spectra used", match.group(1))
                            spectrum.add_attribute(
                                "MS:1003069|number of replicate spectra available", match.group(1))
                        else:
                            spectrum.add_attribute(
                                "ERROR", f"Unable to parse {attributes[attribute]} in {attribute} at E2455")
                            unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)
            #### Expand the Fullname attribute
            elif attribute == "Organism":
                if attributes[attribute] is not None:
                    value = attributes[attribute]
                    value = value.strip('"')

                    if value in species_map:
                        group_identifier = next(group_counter)
                        for item in species_map[value]:
                            spectrum.add_attribute(
                                item[0], item[1], group_identifier)

                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Unable to parse {attributes[attribute]} in {attribute} at E2355")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Otherwise add this term to the list of attributes that we don't know how to handle
            else:
                unknown_terms.append(attribute)

        if "MS:1000888|unmodified peptide sequence" not in spectrum.attribute_dict:
            if "MS:1003061|spectrum name" in spectrum.attribute_dict:
                lookup = spectrum.attribute_dict["MS:1003061|spectrum name"]
                name = spectrum.attributes[lookup["indexes"][0]][1]
                match = re.match("(.+)/(\d+)", name)
                if match:
                    spectrum.add_attribute(
                        "MS:1000888|unmodified peptide sequence", match.group(1))
                    spectrum.add_attribute(
                        "MS:1000041|charge state", match.group(2))

        #### Handle the uninterpretable terms
        for attribute in unknown_terms:
            if attributes[attribute] is None:
                spectrum.add_attribute(
                    "MS:1009900|other attribute name", attribute)
            else:
                group_identifier = next(group_counter)
                spectrum.add_attribute(
                    "MS:1009900|other attribute name", attribute, group_identifier)
                spectrum.add_attribute("MS:1009902|other attribute value",
                                   attributes[attribute], group_identifier)
        return spectrum
