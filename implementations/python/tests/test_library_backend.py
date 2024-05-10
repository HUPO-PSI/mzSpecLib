import math
import os
import unittest

from mzlib.spectrum import Spectrum
from mzlib.backends import (MSPSpectralLibrary, TextSpectralLibrary, JSONSpectralLibrary, SpectronautTSVSpectralLibrary, DIANNTSVSpectralLibrary)
from mzlib.analyte import ANALYTE_MIXTURE_TERM

from .common import datafile


class LibraryBehaviorBase(object):

    def _open_library(self, file_path=None):
        if file_path is None:
            file_path = self.test_file
        lib = self.library_cls(file_path)
        return lib

    def test_sequence_behavior(self):
        lib = self._open_library()
        assert len(lib) == 7
        spec = lib[3]
        assert spec.get_attribute(
            "MS:1003061|library spectrum name") == "AAAAGSTSVKPIFSR/2_0_44eV"

    # TODO: Fix clipping in _buffer_from_stream first
    def test_iteration(self):
        lib = self._open_library()
        assert list(lib) == list(lib.read())


class TestMSPLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.msp")
    library_cls = MSPSpectralLibrary


class MzSpecLibLibraryBehaviorBase(LibraryBehaviorBase):
    def test_interpretation_attribute(self):
        try:
            path = self.test_interpretation_file
        except AttributeError:
            return
        lib = self._open_library(path)
        spec: Spectrum = lib[0]
        interp = spec.interpretations[1]
        assert len(interp.attributes) == 2
        assert interp.get_attribute('MS:1002357|PSM-level probability') == 0.974
        assert interp.get_attribute(ANALYTE_MIXTURE_TERM) == [1, 2]
        mem = interp.get_member_interpretation(1)
        assert len(mem) == 1
        assert mem.get_attribute('MS:XXXXXXX|explained intensity fraction') == 0.287


class TestTextLibrary(unittest.TestCase, MzSpecLibLibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
    library_cls = TextSpectralLibrary
    test_interpretation_file = datafile("complex_interpretations_with_members.mzlb.txt")


class TestJSONLibrary(unittest.TestCase, MzSpecLibLibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.json")
    library_cls = JSONSpectralLibrary
    test_interpretation_file = datafile("complex_interpretations_with_members.mzlb.json")


class TestSpectronautLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("human_serum.head.spectronaut.tsv")
    library_cls = SpectronautTSVSpectralLibrary

    def test_sequence_behavior(self):
        lib = self._open_library()
        assert len(lib) == 9

        spec: Spectrum = lib[0]
        assert spec.name == 'AQIPILR/2'
        assert math.isclose(spec.precursor_mz, 405.7634379)
        assert spec.precursor_charge == 2

        spec = lib[5]
        assert spec.name == 'QELSEAEQATR/2'
        assert math.isclose(spec.precursor_mz, 631.3045807)
        assert spec.get_analyte(1).proteins[0].name == 'CO3_HUMAN'


class TestDIANNTSVLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("phl004_canonical_sall_pv_plasma.head.diann.tsv")
    library_cls = DIANNTSVSpectralLibrary

    def test_sequence_behavior(self):
        lib = self._open_library()
        assert len(lib) == 9

        spec: Spectrum = lib[0]
        analyte = spec.get_analyte(1)
        assert analyte.peptide == 'AAAAAAAAAAAAAAAASAGGK'
        assert spec.name == 'AAAAAAAAAAAAAAAASAGGK2'

