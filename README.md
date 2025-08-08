# mzSpecLib

[![License](https://flat.badgen.net/github/license/HUPO-PSI/mzSpecLib)](https://github.com/HUPO-PSI/SpectralLibraryFormat/blob/master/LICENSE)
[![Open Issues](https://flat.badgen.net/github/open-issues/HUPO-PSI/mzSpecLib)](https://github.com/HUPO-PSI/SpectralLibraryFormat/issues)
[![Open PRs](https://flat.badgen.net/github/open-prs/HUPO-PSI/mzSpecLib)](https://github.com/HUPO-PSI/SpectralLibraryFormat/pulls)
![Contributors](https://flat.badgen.net/github/contributors/HUPO-PSI/mzSpecLib)
![Watchers](https://flat.badgen.net/github/watchers/HUPO-PSI/mzSpecLib)
![Stars](https://flat.badgen.net/github/stars/HUPO-PSI/mzSpecLib)

**HUPO-PSI standardized spectral library format**
mzSpecLib is a formal standard and file format in development at
[HUPO-PSI](http://www.psidev.info/) to store and distribute
spectral libraries/archives. The target main target audience for this format are
the developers of spectral library search tools and resources.

---

- [mzSpecLib](#mzspeclib)
  - [Introduction](#introduction)
  - [Development](#development)
  - [Reference implementation](#reference-implementation)
    - [Demo Viewer](#demo-viewer)
  - [Guidelines](#guidelines)
  - [Contributing](#contributing)

---

## Introduction

Over past years several file formats have been created to store and disseminate
spectral libraries, such as MSP, X!Hunter binary MGF, BiblioSpec SQLite and
SSL/MS2, SpectraST SPLIB/SPTXT,
[MassBank formats](https://github.com/HUPO-PSI/SpectralLibraryFormat/blob/master/legacy-formats/MassBank.md),
Spectronaut CSV...

Each spectral library provider uses one of these formats. For example,
PeptideAtlas uses splib, PRIDE uses MSP, GPMDB  uses X!Hunter binary MGF, and NIST
uses MSP. Some spectral library search engines support multiple formats, and some
do not, making it difficult to share libraries and compare spectral library
searching tools. In the proteomics community, there has been a long-standing
effort to standardize raw mass spectrometric data and the results of data
analysis, primarily identification. But spectral libraries straddle the boundary
between the two and cannot be adequately served by either efforts.

As there remains much fluidity and disagreement in what information should go
into a spectral library, the format must be flexible enough to fit all the
potential use cases of spectral libraries, and yet retain sufficient structure
for it to be a practically useful standard.


## Development

The new format is being developed by the PSI-MS working group. This repository
is used as a central point of information regarding the format's development:
- The issue tracker is used for discussions and decisions regarding the format
- The repository files contain information about current spectral libraries
  formats, the spectral library controlled vocabulary, the spectral library
  specifications, examples and tools that export/validate and visualize those
  files.
- The mzSpecLib specification is nearly complete and has been resubmitted to the PSI Document Process for final community review. It is hoped to be complete at the end of 2024.

If you do experience problems with the files or have suggestions or comments, please open an [issue](https://github.com/HUPO-PSI/mzSpecLib/issues).

- [mzSpecLib main specification document](https://github.com/HUPO-PSI/mzSpecLib/blob/master/specification/mzSpecLib_specification_v1.0_draft09.docx)
- [mzPAF peak interetation format](https://psidev.info/mzPAF)


## Reference implementation

A reference implementation of the mzSpecLib format is available in the form of a Python package.
Check out the [mzspeclib-py](https://github.com/HUPO-PSI/mzspeclib-py) repository or the
[Python package documentation](https://mzspeclib.readthedocs.io/) for more information.

### Demo Viewer

We have a demo web application for viewing libraries: https://mzspeclib-app-viewer-demo.streamlit.app/

## Guidelines

The mzSpecLib format is flexible, but we have examples showing how to use it in certain ways in

## Contributing

All community input is welcome! Feel free to join the discussions in the [Issue
tracker](https://github.com/HUPO-PSI/SpectralLibraryFormat/issues) or to open a
new issue if you have questions, recommendations or requests. Additionally,
everyone is allowed to post comments in the Google documents or to request full
write access to fully contribute to the specification.
