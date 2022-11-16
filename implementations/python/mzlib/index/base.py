import warnings

from typing import Collection, Union, Any, List


class IndexRecordBase:
    __slots__ = ()

    number: int
    offset: int
    name: str


class IndexBase(Collection):

    @classmethod
    def from_filename(cls, filename, library=None):
        raise NotImplementedError()

    def offset_for(self, record_label) -> int:
        record = self.record_for(record_label)
        return record.offset

    def record_for(self, record_label: Union[int, str]) -> IndexRecordBase:
        record = self.search(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def search(self, i: Union[str, int, slice], **kwargs) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        raise NotImplementedError()

    def add(self, number: int, offset: int, name: str, analyte: Any, attributes=None):
        raise NotImplementedError()

    def commit(self):
        raise NotImplementedError()

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, i: Union[int, str, slice]):
        return self.search(i)

    def __len__(self):
        raise NotImplementedError()

    def __contains__(self):
        raise NotImplementedError()

    def check_names_unique(self) -> bool:
        '''Checks that all indexed spectra have unique
        ``spectrum name`` parameters.

        Returns
        -------
        bool:
            Whether the spectrum names in the index are unique.
        '''
        seen = set()
        for record in self:
            if record.name in seen:
                return False
            seen.add(record.name)
        return True
