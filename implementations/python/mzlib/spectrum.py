from __future__ import print_function

import re
import json

from mzlib.attributes import AttributeManager

#A class that holds data for each spectrum that is read from the SpectralLibrary class


class Spectrum(AttributeManager):

    #### Constructor
    def __init__(self, attributes=None):
        """
        __init__ - SpectrumLibrary constructor

        Parameters
        ----------
        attributes: list
            A list of attribute [key, value (, group)] sets to initialize to.

        """
        super(Spectrum, self).__init__(attributes)
        self.peak_list = []

    def write(self, format="text"):
        """
        write - Write out the spectrum in any of the supported formats
        """

        #### Set a buffer to fill with string data
        buffer = ''

        #### Make the format string lower case to facilitate comparisons
        format = format.lower()

        #### If the format is text
        if format == "text":
            from mzlib.backends.text import format_spectrum
            return format_spectrum(self)

        #### If the format is TSV
        elif format == "tsv" or format == "csv":

            #### Set the appropriate delimiter for format
            delimiter = "\t"
            if format == "csv":
                delimiter = ","

            #### Create the header line and columns line
            buffer += "# Spectrum attributes\n"
            buffer += delimiter.join("cv_param_group", "accession", "name", "value_accession", "value\n")

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
                    raise ValueError(
                        f"Unsupported number of items in attribute: {attribute}")
                components = key.split('|',1)
                if len(components) == 2:
                    accession,name = components
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    raise ValueError(f"Unsupported number of items in components: {components}")
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
                    raise ValueError(
                        f"Unsupported number of items in components: {components}")

                #### Create the data line
                buffer += delimiter.join(map(str, [cv_param_group,accession,name,value_accession,value]))+"\n"

            #### Create the header line and columns line
            buffer += "# Peak list\n"
            buffer += "mz\tintensity\tinterpretation\n"

            #### Write out the peak list
            for peak in self.peak_list:
                mz,intensity,interpretation = peak
                if format == "csv" and ',' in str(interpretation):
                    interpretation = '"' + interpretation + '"'
                buffer += delimiter.join(map(str, [mz,intensity,interpretation]))+"\n"

        #### If the format is JSON
        elif format == "json":
            mzs = []
            intensities = []
            interpretations = []
            for peak in self.peak_list:
                mzs.append(peak[0])
                intensities.append(peak[1])
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
                    raise ValueError(
                        f"Unsupported number of items in attribute: {attribute}")
                components = key.split('|',1)
                if len(components) == 2:
                    accession,name = components
                    reformed_attribute['accession'] = accession
                    reformed_attribute['name'] = name
                else:
                    print("ERROR: Unsupported number of items in components")
                    print(components)
                    raise ValueError(
                        f"Unsupported number of items in components: {components}")
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
                    raise ValueError(
                        f"Unsupported number of items in components: {components}")
                attributes.append(reformed_attribute)

            spectrum = { "attributes": attributes, "mzs": mzs, "intensities": intensities,
                "interpretations": interpretations }
            buffer = json.dumps(spectrum,sort_keys=True,indent=2)

        #### Otherwise we don't know this format
        else:
            raise ValueError(f"ERROR: Unrecogized format '{format}'")

        return(buffer)
