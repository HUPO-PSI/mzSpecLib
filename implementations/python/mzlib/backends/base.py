from ..spectrum import Spectrum


class SpectralLibraryBackendBase(object):

    def __init__(self, filename):
        self.filename = filename

    def _new_spectrum(self):
        return Spectrum()

    def get_spectrum(self, spectrum_number=None, spectrum_name=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys):
        raise NotImplementedError()