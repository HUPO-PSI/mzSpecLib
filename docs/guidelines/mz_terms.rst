Specifying ion m/z
------------------

.. csv-table:: M/Z Controlled Vocabulary Terms
   :file: mz_terms.csv
   :widths: 30, 70
   :header-rows: 1
   :name: m/z terms

Selected Ion M/Z
================

The :title-reference:`selected ion m/z` term is intended to be used in :title-reference:`mzSpecLib` in the same manner as in :title-reference:`mzML`. It refers to the m/z value of the precursor that the mass spectrometer targets for isolation and fragmentation, and as such, reflects the assumption that this value is not refined by some other method that the implementation used, and carries with it the same uncertainty. Although in typical mzML files, the value given for :title-reference:`selected ion m/z` is the monoisotopic m/z of the precursor, the definition of the term does not specify it as such. Therefore, for clarity, the term :title-reference:`experimentally determined monoisotopic m/z` is recommended to be used in :title-reference:`mzSpecLib` to refer to the precursor monoisotopic m/z value unambiguously. 

Because it is not directly related to an `Analyte`'s identity, this attribute should be listed in the `Spectrum` section.

Selected Precursor M/Z
======================

The :title-reference:`selected precursor m/z` term should not be used in :title-reference:`mzSpecLib` due to the ambiguity of its definition.

Experimental Precursor Monoisotopic M/Z
=======================================

The :title-reference:`experimental precursor monoisotopic m/z` term was created for :title-reference:`mzSpecLib`. Unlike :title-reference:`selected ion m/z`, which is usually the value copied from the same field in the mzML file, this m/z value might be refined in some way by the implementation to correspond to the monoisotopic m/z of the analyte being fragmented, although it may still deviate from the theoretical monoisotopic m/z due to the limited precision of m/z measurements. This term is the recommended term to encode the experimentally determined precursor monoisotopic m/z of the library spectra in :title-reference:`mzSpecLib`.

Because the determination might be done with information about an `Analyte`, this attribute could appear under either the `Spectrum` (which is more typical) or the `Interpretation` or `InterpretationMember` sections.


Isolation Window Target M/Z
===========================

The :title-reference:`isolation window target m/z` term refers to the reference m/z value around which the mass spectrometer targets for fragmentation to produce the library spectrum. In data-dependent acquisition, it is usually the location of an intense precursor ion peak, which is not necessarily the monoisotopic ion peak. Together with the terms :title-reference:`isolation window upper offset` and :title-reference:`isolation window lower offset`, it can be used to describe the isolation window.  

This attribute should be listed in the `Spectrum` section.


Theroetical Monoisotopic M/Z
=============================

The :title-reference:`theoretical monoisotopic m/z` term refers to the monoisotopic m/z calculated from dividing the theoretical monoisotopic mass of adduct ion (i.e. including the charge carrier(s)) by the charge state. As such, it is completely determined by the identity of the analyte. Therefore, this term should be used exclusively in the `Analyte` section.  

Theroetical Average M/Z
=============================

The :title-reference:`theoretical average m/z` term refers to the average m/z calculated from dividing the theoretical average mass of adduct ion (i.e. including the charge carrier(s)) by the charge state. The average mass refers to the weighted average based on the natural abundance of the isotopes of atoms that make up the analyte. As such, it is completely determined by the identity of the analyte. Therefore, this term should be used exclusively in the `Analyte` section.  




