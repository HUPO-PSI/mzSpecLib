import unittest

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.validate import validator

from .common import datafile


class TestLibrarySemanticValidator(unittest.TestCase):
    def open_library(self):
        test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
        library = SpectrumLibrary(filename=test_file)
        return library

    def validate_library(self, valid):
        library = self.open_library()
        return valid.validate_library(library)

    def test_validate_base(self):
        valid = validator.get_validator_for("base")
        assert self.validate_library(valid)

    def test_validate_peptide(self):
        valid = validator.get_validator_for("base")
        valid = valid.chain(validator.get_validator_for("peptide"))
        assert self.validate_library(valid)

    def test_validate_silver(self):
        valid = validator.get_validator_for("base")
        valid = valid.chain(validator.get_validator_for("silver"))
        assert not self.validate_library(valid)

    def test_validate_peak_annotations(self):
        valid = validator.get_validator_for("base") | validator.get_object_validator_for("peak_annotations")
        assert self.validate_library(valid)

        test_file = datafile("bad_peak_annotations.mzlb.txt")
        library = SpectrumLibrary(filename=test_file)
        assert not valid.validate_library(library)
