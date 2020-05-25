import os
import unittest
import tempfile

from mzlib.backends import MSPSpectralLibrary, TextSpectralLibrary, JSONSpectralLibrary

from .common import datafile


class LibraryBehaviorBase(object):

    def _open_library(self):
        lib = self.library_cls(self.test_file)
        return lib

    def test_sequence_behavior(self):
        lib = self._open_library()
        assert len(lib) == 7
        spec = lib[3]
        assert spec.get_attribute(
            "MS:1003061|spectrum name") == "AAAAGSTSVKPIFSR/2_0_44eV"

    # TODO: Fix clipping in _buffer_from_stream first
    # def test_iteration(self):
    #     lib = self._open_library()
    #     assert list(lib) == list(lib.read())


class TestMSPLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.msp")
    library_cls = MSPSpectralLibrary


class TestTextLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
    library_cls = TextSpectralLibrary


class TestJSONLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.json")
    library_cls = JSONSpectralLibrary
