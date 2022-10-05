from .msp import MSPSpectralLibrary as _MSPSpectralLibrary


class SPTXTSpectralLibrary(_MSPSpectralLibrary):
    file_format = "sptxt"
    format_name = "sptxt"


