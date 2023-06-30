import warnings

from typing import Collection, Iterator, Optional, Union, Any, List


class IndexRecordBase:
    __slots__ = ()

    number: int
    index: int
    offset: int
    name: str


class IndexBase(Collection):

    @classmethod
    def from_filename(cls, filename, library=None):
        raise NotImplementedError()

    def offset_for(self, record_label) -> int:
        record = self.record_for(record_label)
        return record.offset

    def offset_for_cluster(self, record_label) -> int:
        record = self.record_for_cluster(record_label)
        return record.offset

    def record_for(self, record_label: Union[int, str]) -> IndexRecordBase:
        record = self.search(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def record_for_cluster(self, record_label: int) -> IndexRecordBase:
        record = self.search_clusters(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def search(self, i: Union[str, int, slice], **kwargs) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        raise NotImplementedError()

    def search_clusters(self, i: Optional[Union[int, slice]]=None, **kwargs) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        raise NotImplementedError()

    def add(self, number: int, offset: int, name: str, analyte: Any, attributes=None):
        """
        Add a new entry to the spectrum index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        name : str,
            A text identifier for this spectrum.
        analyte : str, optional
            A text representation of the analyte for that record
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        raise NotImplementedError()

    def add_cluster(self, number: int, offset: int, attributes=None):
        """
        Add a new entry to the spectrum index.

        Parameters
        ----------
        number : int
            A numerical identifier for this spectrum.
        offset : int
            The offset in the file to reach the spectrum (in bytes if appropriate)
        attributes : Dict[str, Any], optional
            A key-value pair collection of this record, currently not supported.
        """
        raise NotImplementedError()

    def commit(self):
        """
        Commit any index state to disk, if this index supports persistence.

        Has no effect on index types that do not have a persistence functionality.
        """
        raise NotImplementedError()

    def iter_clusters(self) -> Iterator[IndexRecordBase]:
        raise NotImplementedError()

    def iter_spectra(self) -> Iterator[IndexRecordBase]:
        for i in range(len(self)):
            yield self[i]

    def _get_by_index(self, i: Union[int, slice]) -> Union[IndexRecordBase, List[IndexRecordBase]]:
        raise NotImplementedError()

    def __iter__(self):
        return self.iter_spectra()

    def __getitem__(self, i: Union[int, str, slice]):
        return self.search(i)

    def __len__(self):
        raise NotImplementedError()

    def __contains__(self, key) -> bool:
        try:
            hit = self.search(key)
            return True
        except (KeyError, IndexError, ValueError):
            return False

    def check_names_unique(self) -> bool:
        """
        Checks that all indexed spectra have unique
        ``spectrum name`` parameters.

        Returns
        -------
        bool:
            Whether the spectrum names in the index are unique.
        """
        seen = set()
        for record in self:
            if record.name in seen:
                return False
            seen.add(record.name)
        return True
