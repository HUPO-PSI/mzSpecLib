#!/usr/bin/env python3
# pragma: no cover
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re

class UniversalSpectrumIdentifier(object):

    # usi object takes usiStr an automatically parses it and stores attributes
    # usi objects can still exist even if the usi str is incorrect.
    # it will simply show where the error in the string is

    def __init__(self, usi=None):

        self.is_valid = False
        self.usi = usi
        self.dataset_identifier = None
        self.dataset_subfolder = None
        self.ms_run_name = None
        self.indexFlag = None
        self.index = None
        self.interpretation = None
        self.peptidoform = None
        self.charge = None
        self.provenance_identifier = None
        self.error = 0
        self.error_code = None
        self.error_message = None
        self.warning_message = None

        if usi:
            self.parse(verbose=None)


    # Attributes:
    #   usi
    #   dataset_identifier
    #   dataset_subfolder
    #   ms_run_name
    #   indexFlag
    #   index
    #   interpretation
    #   peptidoform
    #   charge


    #### Set the error state with supplied information
    def set_error(self,error_code,error_message):
        self.error_code = error_code
        self.error_message = error_message
        self.is_valid = False
        #print(f"ERROR: {error_code}: {error_message}")

    # parses USI string
    def parse(self, verbose):
        verboseprint = print if verbose else lambda *a, **k: None
        verboseprint("\nINFO: Parsing USI string '" + self.usi + "'")
        elementOffset = 0
        offset = 0
        if self.usi.startswith("mzspec:"):
            self.usiMzspec = self.usi[len("mzspec:"):]
        else:
            self.set_error("MissingPrefix","USI does not begin with prefix 'mszpec:'")
            return(self)

        # creates list of potential usi attributes
        elements = self.usiMzspec.split(":")
        nElements = len(elements)
        # print(elements)
        # print(nElements)

        # checks if usi has at least 4 colon-separated fields
        if nElements < 4:
            self.set_error("InsufficientComponents","USI string does not have the minimum required 4 colon-separated fields after mzspec:")
            return(self)

        offset = elementOffset

        # dataset_identifier field
        self.dataset_identifier = elements[offset]
        if self.dataset_identifier is None:
            self.set_error("EmptyDatasetIdentifier","Dataset identifier is empty. Not permitted.")
        elif self.dataset_identifier.startswith("PXD") or self.dataset_identifier.startswith("PXL"):
            self.dataset_identifier = elements[offset]
        else:
            self.set_error("UnrecognizedDatasetIdentifier","Dataset identifier unknown. Not permitted.")
            self.error += 1

        elementOffset += 1
        offset = elementOffset
        nextField = elements[offset]
        offsetShift = 0

        #### The old style thing now clashes with the interest in leaving the version out of PXLs
        # empty dataset_subfolder
        #if nextField == '':
        #    self.warning_message = "old style. empty is ok. Empty dataset_subfolder probably."
        #    offsetShift = 1

        offset = elementOffset + offsetShift
        self.ms_run_name = elements[offset]

        if self.ms_run_name:
            verboseprint("MS run equals " + self.ms_run_name)
        else:
            if self.dataset_identifier.startswith("PXL"):
                verboseprint("For PXL identifiers, the version may be blank, okay.")
            else:
                self.set_error("MissingMSRunIdentifier","MS Run identifier empty. Not permitted.")
                self.error += 1

        elementOffset += 1
        offset = elementOffset + offsetShift
        self.indexFlag = elements[offset]
        # print("check " + self.indexFlag)
        # does indexFlag exist?
        if self.indexFlag:
            # is it scan or index
            if self.indexFlag == "scan" or self.indexFlag == "index":
                verboseprint("indexFlag is OK.")
            # is there potentially some weird colon escaping in the msRun name?
            else:
                potentialOffsetShift = offsetShift
                appendStr = ""
                repaired = False

                # fix colon escaping if it exists
                while elementOffset + potentialOffsetShift < nElements:
                    # go until program finds 'scan' or 'index' index flag types
                    if elements[elementOffset + potentialOffsetShift].startswith("scan") or elements[
                        elementOffset + potentialOffsetShift].startswith("index"):
                        self.indexFlag = elements[elementOffset + potentialOffsetShift]
                        self.ms_run_name += appendStr
                        offsetShift = potentialOffsetShift
                        repaired = True
                        break
                    appendStr += ":" + elements[elementOffset + potentialOffsetShift]

                    potentialOffsetShift += 1

                # colon escape fixed and msRun field updated
                if repaired:
                    verboseprint("Unescaped colon in msRun name. Hopefully taken care of. Please fix this")
                    verboseprint("msRun name revised to '{}'".format(self.ms_run_name))

                # no 'scan' or 'index' fields found later. assume broken index flag
                else:
                    self.error += 1
                    self.set_error("InvalidIndexType","Index type invalid. Must be 'scan' or 'index'")

        # no index flag
        else:
            self.error += 1
            self.set_error("MissingIndexType","Index flag empty! Not permitted.")

        elementOffset += 1
        offset = offsetShift + elementOffset

        # index for index flag if flag is valid. useless if index flag is invalid
        self.index = elements[offset]
        if self.index:
            verboseprint("Index is " + self.index)
        else:
            self.set_error("MissingIndex","Index number empty! Not permitted.")
            self.error += 1

        elementOffset += 1
        offset = elementOffset + offsetShift

        # if statement check to see if the USI even has an interpretation field
        if offset < nElements:
            self.interpretation = elements[offset]
            self.peptidoform = ''
            self.charge = ''
            if self.interpretation and self.interpretation != '':
                find = re.match(r"^\s*(.+)\/(\d+)\s*$", self.interpretation)
                # match
                if find:
                    # subfields of interpretation
                    self.peptidoform = find.group(1)
                    self.charge = find.group(2)
                    verboseprint("Interpreted peptidoform = {}, charge = {}".format(self.peptidoform, self.charge))
                else:
                    self.set_error("MissingIndex","Unable to parse interpretation {} as peptidoform/charge".format(self.interpretation))

            else:
                verboseprint("Interpretation field not provided. OK.")

        # provenance identifier
        if offset < nElements:
            self.provenance_identifier = elements[offset]
            verboseprint("Provenance Identifier = ".format(self.provenance_identifier))
        # returns count of errors found in usi. useful for checking if the entire identifier is valid.

        if self.error == 0:
           self.is_valid = True
        # errors found in usi
        else:
            verboseprint("Number of errors: " + str(self.error))
            self.is_valid = False
            verboseprint("ERROR: Invalid USI " + self.usi)

    # prints out USI attributes
    def show(self):
        print("USI: " + str(self.usi))
        print("is_valid: " + str(self.is_valid))
        print("Dataset Identifier: " + str(self.dataset_identifier))
        print("Dataset Subfolder: " + str(self.dataset_subfolder))
        print("MS run name: " + str(self.ms_run_name))
        print("Index flag: " + str(self.indexFlag))
        print("Index: " + str(self.index))
        print("Peptidoform: " + str(self.peptidoform))
        print("Charge: " + str(self.charge))


# If this class is run from the command line, perform a short little test to see if it is working correctly
def run_tests():
    testUSIs = [
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951"],
        ["invalid", "PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951"],
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2"],
        ["invalid", "mzspec:PASS002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2"],
        ["invalid", "mzspec"],
        ["invalid", "mzspec:"],
        ["invalid", "mzspec:PXD001234"],
        ["invalid", "mzspec:PXD001234:00261_A06_P001564_B00E_A00_R1:scan"],
        ["valid", "mzspec:PXD001234:00261_A06_P001564_B00E_A00_R1:index:10951"],
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[Phospho]IDELVISK/2"],
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951:PEPT[+79]IDELVISK/2"],
        ["valid", "mzspec:PXD001234:Dilution1:4:scan:10951"],
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:test1:scan:10951:PEPT[Phospho]IDELVISK/2"],
        ["valid", "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1\\:test1:scan:10951:PEPT[Phospho]IDELVISK/2"],
    ]
    testUSIsValid = []

    # Loop over each test USI, parse it, and determine if it is valid or not, and print the index number
    print("Testing example USIs:")
    for usiSet in testUSIs:
        expectedStatus = usiSet[0]
        usiStr = usiSet[1]

        # Create a new UniversalSpectrumIdentifier object
        # made the USI object itself take a string so that parse does not need to be called explicitly
        usi = UniversalSpectrumIdentifier(usiStr)
        expected_validity = True
        if expectedStatus == 'invalid':
            expected_validity = False
        else:
            expectedStatus = 'valid  '

        status = 'PASS'
        if usi.is_valid is not expected_validity:
            status = 'FAIL'

        response = usi.is_valid
        testUSIsValid.append(response)
        print(f"{status}\texpected {expectedStatus}\t{usiStr}")
    # check to see if parsing is correct
    #print(testUSIsValid)


#### A very simple example of using this class
def example():
    usiStr = "mzspec:PXD002437:00261_A06_P001564_B00E_A00_R1:scan:10951"
    usiStr = "mzspec:PXL000001:01-09-2015:index:500"
    usiStr = "mzspec:PXL000001::index:500"
    usi = UniversalSpectrumIdentifier(usiStr)
    #usi.parse(verbose=1)
    usi.show()


#### If class is run directly
def main():
    #example()
    run_tests()


if __name__ == "__main__": main()
