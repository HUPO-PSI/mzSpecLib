try:
    from collections.abc import (MutableMapping, Mapping)
except ImportError:
    from collections import (MutableMapping, Mapping)

import textwrap

from mzlib.attributes import AttributedEntity, AttributeManager


FIRST_ANALYTE_KEY = '1'
FIRST_INTERPRETATION_KEY = '1'


class _AnalyteMappingProxy(Mapping):
    def __init__(self, parent):
        self.parent = parent

    def __getitem__(self, key):
        for group_id, group in self.parent.items():
            try:
                return group[key]
            except KeyError:
                pass
        raise KeyError(key)

    def __iter__(self):
        keys = set()
        for group_id, group in self.parent.items():
            k_ = set(group.keys())
            assert not (keys & k_)
            keys.update(k_)
        return iter(keys)

    def __len__(self):
        n = 0
        for group_id, group in self.parent.items():
            n += len(group)
        return n

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class InterpretationCollection(MutableMapping):
    __slots__ = ('interpretations', )

    def __init__(self, interpretations=None):
        if interpretations is None:
            interpretations = {}
        self.interpretations = interpretations

    def get_interpretation(self, interpretation_id) -> 'Interpretation':
        return self.interpretations[str(interpretation_id)]

    def add_interpretation(self, interpretation: 'Interpretation'):
        self.set_interpretation(str(interpretation.id), interpretation)

    def set_interpretation(self, key, interpretation: 'Interpretation'):
        self.interpretations[str(key)] = interpretation

    def __getitem__(self, key) -> 'Interpretation':
        return self.get_interpretation(key)

    def __setitem__(self, key, value: 'Interpretation'):
        self.set_interpretation(key, value)

    def __delitem__(self, key):
        del self.interpretations[str(key)]

    def __len__(self):
        return len(self.interpretations)

    def __contains__(self, key):
        return key in self.interpretations

    def __iter__(self):
        return iter(self.interpretations)

    def keys(self):
        return self.interpretations.keys()

    @property
    def analytes(self) -> '_AnalyteMappingProxy':
        return _AnalyteMappingProxy(self)

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"



class Interpretation(AttributedEntity, MutableMapping):
    __slots__ = ('id', 'analytes', )

    id: str
    analytes: dict

    def __init__(self, id, attributes=None, analytes=None):
        self.id = str(id)
        self.analytes = analytes or {}
        super(Interpretation, self).__init__(attributes)

    def get_analyte(self, analyte_id) -> 'Analyte':
        return self.analytes[str(analyte_id)]

    def add_analyte(self, analyte: 'Analyte'):
        self.set_analyte(analyte.id, analyte)

    def set_analyte(self, key, analyte: 'Analyte'):
        self.analytes[str(key)] = analyte

    def remove_analyte(self, analyte_id):
        del self.analytes[str(analyte_id)]

    def __getitem__(self, key) -> 'Analyte':
        return self.get_analyte(key)

    def __setitem__(self, key, value: 'Analyte'):
        self.set_analyte(key, value)

    def __delitem__(self, key):
        self.remove_analyte(key)

    def __iter__(self):
        return iter(self.analytes)

    def __len__(self):
        return len(self.analytes)

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class Analyte(AttributeManager):
    __slots__ = ('id', )

    id: str

    def __init__(self, id, attributes=None):
        self.id = str(id)
        super(Analyte, self).__init__(attributes)

    def __repr__(self):
        template = f"{self.__class__.__name__}(id={self.id}, "
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[])"
            return template
        template += "[\n%s])" % textwrap.indent(',\n'.join(lines), ' ' * 2)
        return template
