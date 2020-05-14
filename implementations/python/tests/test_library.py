import os
import unittest
import tempfile

from mzlib.spectrum_library import SpectrumLibrary

from .common import datafile
from .test_library_backend import LibraryBehaviorBase


class TestSpectrumLibrary(unittest.TestCase, LibraryBehaviorBase):
    library_cls = SpectrumLibrary
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")

    def _open_library(self):
        lib = self.library_cls(filename=self.test_file)
        return lib
