import os

from mzlib import SpectrumLibrary

workspace = os.path.dirname(__file__)

source = SpectrumLibrary(filename=os.path.join(
    workspace, "chinese_hamster_hcd_selected_head.msp"))
print(source)
annots = set()
for spec in source:
    print(spec.name)
    for peak in spec.peak_list:
        for a in peak[2]:
            annots.add(a)

with open(os.path.join(workspace, "annotations.txt"), 'wt') as fh:
    fh.write("\n".join(sorted(map(str, annots))))
