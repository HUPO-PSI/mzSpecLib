import re

annotation_pattern = re.compile(r"""
(?:(?:(?P<series>[bycz]\.?)(?P<ordinal>\d+))|
   (?P<series_internal>[m](?P<internal_start>\d+):(?P<internal_end>\d+))|
   (?P<precursor>p)|
   (?P<immonium>[ARNDCEQGHKMFPSTWYV]|IL|LI|I|L)|
   (?P<reporter>r(?P<reporter_mass>\d+(?:\.\d+)))|
   (?:_(?P<external_ion>[^\s,/]+))
)
(?P<neutral_loss>[+-](?:[A-Za-z0-9]+))?
(?:(?P<isotope>[+-]\d*i))?
(?:\^(?P<charge>[+-]?\d+))?
(?:@(?P<analyte_reference>[^/\s]+))?
""", re.X)

# At the time of first writing, this pattern could be translated into the equivalent
# ECMAScript compliant regex:
# (?:(?:(?<series>[bycz])(?<ordinal>\d+))|(?<series_internal>[m](?<internal_start>\d+):(?<internal_end>\d+))|(?<precursor>p)|(?<immonium>[ARNDCEQGHKMFPSTWYV]|IL|LI|I|L)|(?<reporter>r(?<reporter_mass>\d+(?:\.\d+)))|(?:_(?<external_ion>[^\s,/]+)))(?<neutral_loss>[+-](?:[A-Za-z0-9]+))?(?:(?<isotope>[+-]\d*i))?(?:\^(?<charge>[+-]?\d+))?(?:@(?<analyte_reference>[^/\s]+))?
# Line breaks not introduced to preserve syntactic correctness.

