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

- [Introduction](#introduction)
- [Development](#development)
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
- Additionally, Google Docs are used to quickly and efficiently collaborate
  on the specification documents. Links to these documents are listed here:
  - [mzSpecLib main specification document](https://docs.google.com/document/d/1l87lIyKTy2ti5yU7aqsLr7uX5jIU1dO7gEzyqWD2uQA/edit?usp=sharing)
  - [Metadata and CV term overview](https://drive.google.com/file/d/1rN5DJSowp2micxlwJQlPxlv39ZiaLEfv/edit?usp=sharing)
  - [mzSpecLib peak interetation format](https://docs.google.com/document/d/1yEUNG4Ump6vnbMDs4iV4s3XISflmOkRAyqUuutcCG2w/edit?usp=sharing)
  - [mzSpecLib general data model schematic](https://drive.google.com/file/d/1OVh5ATfKXA77pM4CYzRfdupeRGu3vt5c/view?usp=sharing)


Currently, the project's progress can be split up into the development of the
[main specification](https://github.com/HUPO-PSI/mzSpecLib/tree/master/specification)
and into a [Python implementation](https://github.com/HUPO-PSI/mzSpecLib/tree/master/implementations/python)
of this specification.


## Contributing

All community input is welcome! Feel free to join the discussions in the [Issue
tracker](https://github.com/HUPO-PSI/SpectralLibraryFormat/issues) or to open a
new issue if you have questions, recommendations or requests. Additionally, 
everyone is allowed to post comments in the Google documents or to request full
write access to fully contribute to the specification.
