from __future__ import print_function

import re
import textwrap
import json

from mzlib.attributes import AttributeManager

#A class that holds data for each spectrum that is read from the SpectralLibrary class


class Spectrum(AttributeManager):

    #### Constructor
    def __init__(self, attributes=None, peak_list=None, analytes=None):
        """
        __init__ - SpectrumLibrary constructor

        Parameters
        ----------
        attributes: list
            A list of attribute [key, value (, group)] sets to initialize to.
        """
        if peak_list is None:
            peak_list = []
        if analytes is None:
            analytes = []
        super(Spectrum, self).__init__(attributes)
        self.peak_list = peak_list
        self.analytes = analytes

    def __eq__(self, other):
        result = super(Spectrum, self).__eq__(other)
        if result:
            result = self.peak_list == other.peak_list
        if result:
            result = self.analytes == other.analytes
        return result

    def __repr__(self):
        template = f"{self.__class__.__name__}("
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[], "
        else:
            template += "[\n%s], " % textwrap.indent(',\n'.join(lines), ' ' * 2)
        lines = list(map(str, self.peak_list))
        if not lines:
            template += "peak_list=[])"
        else:
            template += "peak_list=[\n%s])" % textwrap.indent(
                ',\n'.join(lines), ' ' * 2)
        return template

    def __str__(self):
        return self.write("text")

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
            from mzlib.backends.json import format_spectrum
            return format_spectrum(self)

        #### Otherwise we don't know this format
        else:
            raise ValueError(f"ERROR: Unrecogized format '{format}'")

        return(buffer)
