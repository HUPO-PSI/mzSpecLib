class IndexBase(object):

    @classmethod
    def from_filename(cls, filename, library=None):
        raise NotImplementedError()

    def offset_for(self, record_label):
        record = self.record_for(record_label)
        return record.offset

    def record_for(self, record_label):
        record = self.search(record_label)
        if isinstance(record, list):
            warnings.warn(
                f"Multiple records found for {record_label}, using the first")
            record = record[0]
        return record

    def search(self, i, **kwargs):
        raise NotImplementedError()

    def add(self, number, offset, name, analyte, attributes=None):
        raise NotImplementedError()

    def commit(self):
        raise NotImplementedError()
