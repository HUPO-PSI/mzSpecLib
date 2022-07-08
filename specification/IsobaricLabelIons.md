# Isobaric Label Peak Names

In mzSpecLib library peak annotations, reporter ions and other peaks related to reporter ions will be labeled with
the following conventions. A useful resource for currently available reporter ions is available at the UWPR:
https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php. Additional information about TMT is
available in the Thermo Fisher TMT brochure:
https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0015866_2162600_TMT_MassTagging_UG.pdf

# TMT peaks

- TMT126
- TMT127N
- TMT127C
- TMT128N
- TMT128C
- TMT129N
- TMT129C
- TMT130N
- TMT130C
- TMT131N
- TMT131C
- TMT0Nterm
- TMT2Nterm
- TMT6Nterm
- TMTproNterm
- TMT132N
- TMT132C
- TMT133N
- TMT133C
- TMT134N

# iTRAQ peaks

- iTRAQ113
- iTRAQ114
- iTRAQ115
- iTRAQ116
- iTRAQ117
- iTRAQ118
- iTRAQ119
- iTRAQ121
- iTRAQ4Nterm
- iTRAQ8Nterm

The "Nterm" elements refer to the entire N terminal label, including both the reporter and balance group. These commonly fall
off together producing a peak at the full label mass plus a proton.

Each of these may also be used as a neutral (or charged actually) loss. It is also common to to see
- p-[TMT6Nterm]

It is also not unusual to see a series of peaks such as:
- p-[iTRAQ114]
- p-[iTRAQ115]
- p-[iTRAQ116]-H2O
- etc.

wherein the reporter ion falls off with a charge (yielding the reporter ion peak) and the full peptide sequence plus the
balance group is seen as a singly charged ion.


