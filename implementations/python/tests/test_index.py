import os
import unittest
import tempfile

from mzlib.backends import TextSpectralLibrary
from mzlib.index import MemoryIndex, SQLIndex

from .common import datafile

class IndexBehaviorBase(object):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
    index_cls = None

    def make_dir(self):
        dirname = tempfile.mkdtemp("_test", "mzlb_")
        return dirname

    def _open(self):
        return TextSpectralLibrary(self.test_file)

    def _make_index(self, library):
        library.index = self.index_cls()
        library.create_index()
        return library.index

    def test_sequence_behavior(self):
        lib = self._open()
        index = self._make_index(lib)
        assert len(index) == 7
        record = index[3]
        assert record.number == 3
        assert record.name == "AAAAGSTSVKPIFSR/2_0_44eV"


class TestMemoryIndex(unittest.TestCase, IndexBehaviorBase):
    index_cls = MemoryIndex


class TestSQLIndex(unittest.TestCase, IndexBehaviorBase):
    index_cls = SQLIndex

    def _make_index(self, library):
        dirname = self.make_dir()
        index_root = os.path.join(dirname, os.path.basename(self.test_file))
        library.index = self.index_cls(index_root)
        library.create_index()
        return library.index
