Specifying Ion Mass
-------------------

.. csv-table:: Mass Controlled Vocabulary Terms
   :file: mass_terms.csv
   :widths: 30, 70
   :header-rows: 1
   :name: ion mass terms


Molecular Mass
==============

The :title-reference:`molecular mass` term is a generic term used to describe the mass of a molecule in Da. While in typical use it refers to the average mass based on the natural isotope abundances, the term is too vague to be used in the context of MS. Its use in :title-reference:`mzSpecLib` is discouraged in favor of the more explicit term :title-reference:`theroetical neutral average mass`. 


Adduct Ion Mass
===============

The :title-reference:`adduct ion mass` term refers to the mass of the adduct ion which fragments to yield that library spectrum. The mass should include both the neutral analyte and the charge-carrying moeity. However, it does not specify whether it is the monoisotopic mass or the average mass. For clarity, :title-reference:`adduct ion monoisotopic mass` and :title-reference:`adduct ion average mass` should be used instead. 

Adduct Ion Monoisotoptic Mass
=============================

The :title-reference:`adduct ion monoisotopic mass` term refers to the *monoisotopic* mass of the adduct ion which fragments to yield that library spectrum. The mass should include both the neutral analyte and the charge-carrying moeity. This value should be calculated from summing the monoisotopic masses of the constituent atoms/ions, and as such, should be listed under the `Analyte` section.


Adduct Ion Average Mass
=======================

The :title-reference:`adduct ion average mass` term refers to the *average* mass of the adduct ion which fragments to yield that library spectrum. The mass should include both the neutral analyte and the charge-carrying moeity. This value should be calculated from summing the average masses of the constituent atoms/ions, and as such, should be listed under the `Analyte` section.

Theoretical Neutral Mass
========================

The :title-reference:`theoretical neutral mass` term refers to the mass of the neutral analyte. It does not include the charge-carrying moeity. However, it does not specify whether it is the monoisotopic mass or the average mass. For clarity, :title-reference:`theoretical neutral monoisotopic mass` and :title-reference:`theoretical neutral average mass` should be used instead. 

Theoretical Neutral Monoisotopic Mass
=====================================

The :title-reference:`theoretical neutral monoisotopic mass` term refers to the *monoisotopic* mass of the neutral analyte. It does not include the charge-carrying moeity. This value should be calculated from summing the monoisotopic masses of the constituent atoms, and as such, should be listed under the `Analyte` section.

Theoretical Neutral Average Mass
================================

The :title-reference:`theoretical neutral average mass` term refers to the *average* mass of the neutral analyte. It does not include the charge-carrying moeity. This value should be calculated from summing the average masses of the constituent atoms, and as such, should be listed under the `Analyte` section.




