import unittest

from mzlib.backends import (MSPSpectralLibrary)
from mzlib.backends import msp

from .common import datafile


class TestMSPAttributeConversion(unittest.TestCase):
    def test_hierarchy(self):
        flat_maps = (msp.analyte_terms, msp.other_terms, msp.interpretation_terms, msp.interpretation_member_terms)
        handler_maps = (msp.msp_analyte_attribute_handler, msp.msp_spectrum_attribute_handler)
        for handler_map in handler_maps:
            for key in handler_map.mapping:
                for flat_map in flat_maps:
                    assert key not in flat_map
