#!/usr/bin/env python
import io
import tarfile
import pprint

from urllib.request import urlopen

from mzlib.backends import msp
from mzlib import annotation


urls = [
    "https://chemdata.nist.gov/download/peptide_library/libraries/human/HCD/2020_05_19/human_hcd_tryp_good.msp.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/human/HCD/2019_02_14/human_hcd_labelfree_phospho_selected_passed.msp.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/cptaclib/2015/cptac2_human_hcd_itraq_phospho_selected.msp.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/human/HCD/2022_03_16/cptac3_tmt_selected_passed_best.msp.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/sigmaups1/ION_TRAP/2011_05_24/2011_05_24_sigmaups1_consensus_final_true_lib.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/proteins/ion_trap/bsa/strict/2011_04_01/2011_04_01_bsa_consensus_final_true_lib.tar.gz",
    "https://chemdata.nist.gov/download/peptide_library/libraries/cptaclib/2015/cptac2_mouse_hcd_itraq_phospho_selected.msp.tar.gz",
]


def open_url(url: str, buffer_size: int=2 ** 20):
    buffer = io.BytesIO(
        urlopen(url).read(buffer_size)
    )

    arc = tarfile.open(
        fileobj=buffer,
        mode='r:gz')

    ti = arc.firstmember
    fh = arc.extractfile(ti)
    return fh


def test_url(url: str, n_spectra_to_read: int=100, track_unknown_attribute_values: bool=True):
    buffer = open_url(url)
    backend = msp.MSPSpectralLibrary(buffer, read_metadata=False, create_index=False)

    if track_unknown_attribute_values:
        backend.unknown_attributes = msp.UnknownKeyValueTracker()
    for i, spec in enumerate(backend.read()):
        for (mz, inten, annots, aggs) in spec.peak_list:
            for annot in annots:
                if isinstance(annot, annotation.InvalidAnnotation):
                    print(f"Failed to parse {annot} with m/z {mz}")
        if i >= n_spectra_to_read:
            break

    pprint.pprint(backend.summarize_parsing_errors())


def main():
    for url in urls:
        print(f"Testing URL {url}")
        test_url(url)


if __name__ == "__main__":
    main()
