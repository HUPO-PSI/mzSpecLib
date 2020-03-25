import itertools

from ..spectrum import Spectrum


class SpectralLibraryBackendBase(object):

    def __init__(self, filename):
        self.filename = filename

    def _make_counter(self):
        return itertools.count(1)

    def _new_spectrum(self):
        return Spectrum()