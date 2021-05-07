import os
import unittest
import json

from mzlib.annotation import parse_annotation, MassError, IonAnnotationBase

from .common import datafile

try:
    import jsonschema
except ImportError:
    jsonschema = None


DEFAULT_SCHEMA_LOCATION = "../../specification/peak-annotation-format/annotation-schema.json"
SCHEMA_LOCATION = os.environ.get("MZSPECLIB_PEAK_ANNOTATION_SCHEMA_LOCATION", DEFAULT_SCHEMA_LOCATION)
try:
    schema_error = None
    with open(SCHEMA_LOCATION, 'rt') as fh:
        schema = json.load(fh)
except OSError as err:
    schema = None
    schema_error = err


skip_reason = ''
if jsonschema is None:
    skip_reason += ("The `jsonschema` library is not installed, cannot test schema compliance. ")

if schema is None:
    skip_reason += ("Could not accessed `peak-annotation-format.json` from development, cannot test schema compliance. ")

class TestAnnotationParser(unittest.TestCase):
    def test_parse_failfast(self):
        assert parse_annotation("?") == []

    def test_parse_annotation_complex(self):
        base = "b14"
        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14

        base += "-H2O-NH3+[Foo]"

        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']

        base += "+2i"

        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2

        base += "^2"
        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2

        base += "[M+NH4]"
        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2
        assert parsed.adducts == ["M", "NH4"]

        base += "/0.5ppm"

        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2
        assert parsed.adducts == ["M", "NH4"]
        assert parsed.mass_error == MassError(0.5, 'ppm')

        base = "2@" + base

        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2
        assert parsed.adducts == ["M", "NH4"]
        assert parsed.mass_error == MassError(0.5, 'ppm')
        assert parsed.analyte_reference == '2'
        assert parsed == base

        base = base + '*0.05'

        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2
        assert parsed.adducts == ["M", "NH4"]
        assert parsed.mass_error == MassError(0.5, 'ppm')
        assert parsed.analyte_reference == '2'
        assert parsed.confidence == 0.05
        assert parsed == base

        base = '[%s]' % base
        parsed = parse_annotation(base)[0]
        assert parsed.series == 'b'
        assert parsed.position == 14
        assert parsed.neutral_loss == ['-H2O', '-NH3', '[Foo]']
        assert parsed.isotope == 2
        assert parsed.charge == 2
        assert parsed.adducts == ["M", "NH4"]
        assert parsed.mass_error == MassError(0.5, 'ppm')
        assert parsed.analyte_reference == '2'
        assert parsed.confidence == 0.05
        assert parsed.is_auxiliary
        assert parsed == base

    @unittest.skipIf(schema is None or jsonschema is None, skip_reason)
    def test_jsonschema_compliance(self):
        with open(datafile('annotations.txt'), 'rt') as fh:
            for line in fh:
                line = line.strip()
                annot = parse_annotation(line)
                if not annot:
                    raise ValueError("Could not parse %r" % annot)
                annot = annot[0]
                data = annot.to_json()
                jsonschema.validate(instance=data, schema=schema)

    def test_from_json(self):
        with open(datafile('annotations.txt'), 'rt') as fh:
            for line in fh:
                line = line.strip()
                annot = parse_annotation(line)
                if not annot:
                    raise ValueError("Could not parse %r" % annot)
                annot = annot[0]
                data = annot.to_json()
                parsed = IonAnnotationBase.from_json(data)
                # JSON schema enforces analyte_reference but it is implicitly '1'
                # without explicitly being set
                annot.analyte_reference = '1'
                self.assertEqual(annot, parsed)
