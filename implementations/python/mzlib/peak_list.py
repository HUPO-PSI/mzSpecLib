import enum

from pprint import pformat
from numbers import Number
from typing import List, Union, Mapping, Sequence, Iterator

import numpy as np

PeakDtype = np.dtype([('mz', np.float64), ('intensity', np.float32),
                      ('annotation', object), ('aggregation', object)])

PeakType = Mapping[str, Union[float, list]]

class ErrorUnit(enum.Enum):
    Da = 'da'
    PPM = 'ppm'


class PeakList(Sequence[PeakType]):
    peaks: np.ndarray

    def __init__(self, peaks):
        if isinstance(peaks, np.ndarray) and peaks.dtype == PeakDtype:
            self.peaks = peaks
        else:
            self.peaks = np.array([tuple(p) for p in peaks], dtype=PeakDtype)

    def __len__(self):
        return len(self.peaks)

    def __getitem__(self, i):
        return self.peaks[i]

    def __iter__(self) -> Iterator[PeakType]:
        return iter(self.peaks)

    def __repr__(self):
        return f"{self.__class__.__name__}({pformat(self.peaks, indent=2)})"

    def __eq__(self, other):
        if other is None:
            return False
        if not isinstance(other, PeakList):
            other = PeakList(other)
        valid = np.allclose(self.peaks['mz'], other.peaks['mz'])
        if not valid:
            return False
        valid = np.allclose(self.peaks['intensity'], other.peaks['intensity'])
        if not valid:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def find(self, mz, error_tolerance=10, error_unit='ppm') -> List[Sequence[PeakType]]:
        if isinstance(mz, Number):
            mz = [mz]
        mzs = self.peaks['mz']
        ii = np.searchsorted(mzs, mz)
        n = len(self)
        error_unit = ErrorUnit(error_unit)

        outs = []
        if error_unit == ErrorUnit.Da:
            for i, mz_i in zip(ii, mz):
                low = i - 1
                while low >= 0:
                    if abs(mzs[low] - mz_i) < error_tolerance:
                        low -= 1
                    else:
                        low += 1
                        break
                high = i
                while high < n:
                    if abs(mzs[high] - mz_i) < error_tolerance:
                        high += 1
                    else:
                        high -= 1
                        break
                outs.extend(self.peaks[slice(low, high + 1)])
        elif error_unit == ErrorUnit.PPM:
            error_tolerance /= 1e6
            for i, mz_i in zip(ii, mz):
                low = i - 1
                while low >= 0:
                    if abs(mzs[low] - mz_i) / mz_i < error_tolerance:
                        low -= 1
                    else:
                        low += 1
                        break
                high = i
                while high < n:
                    if abs(mzs[high] - mz_i) / mz_i < error_tolerance:
                        high += 1
                    else:
                        high -= 1
                        break
                outs.extend(self.peaks[slice(low, high + 1)])
        return outs
