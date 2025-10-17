Specifying Ion Mass
-------------------

.. csv-table:: Mass Controlled Vocabulary Terms
   :file: mass_terms.csv
   :widths: 30, 70
   :header-rows: 1
   :name: ion mass terms


Molecular Mass
==============

The :title-reference:`molecular mass` term is a generic term used to describe the mass of a molecule in Da. It is too vague to be used in the context of MS. Its use in
:title-reference:`mzSpecLib` is discouraged.


Adduct Ion Mass
==============

The :title-reference:`adduct ion mass` term refers to the mass of the adduct ion which fragments to yield that library spectrum. The mass should include both the neutral analyte and the charge-carrying moeity. It is expected that this value is calculated from summing the *monoisotopic* masses of the constituent atoms/ions, and as such, should be listed under the `Analyte` section. 


Theoretical Neutral Mass
========================

The :title-reference:`theoretical neutral mass` term refers to the mass of the neutral analyte. It does not include the charge-carrying moeity. It is expected that this value is calculated from summing the *monoisotopic* masses of the constituent atoms, and as such, should be listed under the `Analyte` section.

