import sys
import matplotlib
matplotlib.use("agg")

from matplotlib import pyplot as plt

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.draw import draw_spectrum


def main(path, spectrum_key):
    lib = SpectrumLibrary(filename=path)
    spec = lib.get_spectrum(spectrum_number=spectrum_key)
    draw_spectrum(spec)
    plt.savefig(f"{path}.{spectrum_key}.annotated.pdf", bbox_inches='tight')


if __name__ == '__main__':
    try:
        path = sys.argv[1]
        index = sys.argv[2]
        main(path, int(index))
        sys.exit(0)
    except (IndexError, TypeError):
        print("USAGE: <prog> <spectrum library path> <spectrum key>")
        print("\tWrites the annotated spectrum to <spectrum library path>.<spectrum key>.annotated.pdf")
        sys.exit(1)
