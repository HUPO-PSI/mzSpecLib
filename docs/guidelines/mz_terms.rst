Specifying ion m/z
------------------

.. csv-table:: M/Z Controlled Vocabulary Terms
   :file: mz_terms.csv
   :widths: 30, 70
   :header-rows: 1
   :name: m/z terms

Selected Ion M/Z
================

The :title-reference:`selected ion m/z` term is intended to be used in
:title-reference:`mzSpecLib` in the same manner as in :title-reference:`mzML`.
It reflects the assumption that this value is not refined by some other
method that the implementation used, and carries with it the same uncertainty.
Because it is not directly related to an `Analyte`'s identity, this attribute
belongs to the `Spectrum` alone.

Experimentally Determined Monoisotopic M/Z
==========================================

The :title-reference:`experimentally determined monoisotopic m/z` term
was created for :title-reference:`mzSpecLib` to indicate that this m/z
value was refined in some way by the implementation. This implies that
the value is expected to be *correct* in-so-far as the implementation
can tell.

Because the determination might be done with information about
an `Analyte`, this attribute could appear under either the `Spectrum` or
the `Interpretation` or `InterpretationMember` sections.
