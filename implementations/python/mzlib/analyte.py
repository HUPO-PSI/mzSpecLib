import warnings

try:
    from collections.abc import (MutableMapping, Mapping)
except ImportError:
    from collections import (MutableMapping, Mapping)

from typing import Iterable, KeysView, ItemsView, Optional, ValuesView, Dict

from pyteomics import proforma

from mzlib.attributes import AttributedEntity, IdentifiedAttributeManager, AttributeManagedProperty, AttributeProxy, AttributeGroupFacet


FIRST_ANALYTE_KEY = '1'
FIRST_INTERPRETATION_KEY = '1'

ANALYTE_MIXTURE_TERM = "MS:1003163|analyte mixture members"
CHARGE_STATE = "MS:1000041|charge state"
PROFORMA_ION = "MS:1003270|proforma peptidoform ion notation"
PROFORMA_SEQ = "MS:1000889|proforma peptidoform sequence"


class _AnalyteMappingProxy(Mapping):
    parent: Mapping

    def __init__(self, parent):
        self.parent = parent

    def __getitem__(self, key):
        for group_id, group in self.parent.items():
            if key in group:
                return group[key]
        raise KeyError(key)

    def __iter__(self):
        keys = set()
        for group_id, group in self.parent.items():
            k_ = set(group.keys())
            keys.update(k_)
        return iter(keys)

    def __len__(self):
        return sum([len(v) for v in self.parent.values()])

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class InterpretationCollection(MutableMapping):
    __slots__ = ('interpretations', )

    interpretations: Dict[str, 'Interpretation']

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

    def keys(self) -> KeysView[str]:
        return self.interpretations.keys()

    def values(self) -> ValuesView['Interpretation']:
        return self.interpretations.values()

    def items(self) -> ItemsView[str, 'Interpretation']:
        return self.interpretations.items()

    @property
    def analytes(self) -> '_AnalyteMappingProxy':
        return _AnalyteMappingProxy(self)

    def __repr__(self):
        d = dict(self)
        return f"{self.__class__.__name__}({d})"


class Interpretation(AttributedEntity, MutableMapping):
    __slots__ = ('id', 'analytes', 'member_interpretations')

    id: str
    analytes: Dict[str, 'Analyte']
    member_interpretations: Dict[str, 'InterpretationMember']

    def __init__(self, id, attributes: Iterable = None, analytes: Dict = None, member_interpretations: Dict = None):
        self.id = str(id)
        self.analytes = analytes or {}
        self.member_interpretations = member_interpretations or {}
        super(Interpretation, self).__init__(attributes)

    def _update_mixture_members_term(self):
        value = sorted(map(int, self.analytes.keys()))
        self.replace_attribute(ANALYTE_MIXTURE_TERM, value)

    def get_analyte(self, analyte_id) -> 'Analyte':
        return self.analytes[str(analyte_id)]

    def add_analyte(self, analyte: 'Analyte'):
        self.set_analyte(analyte.id, analyte)
        self._update_mixture_members_term()

    def set_analyte(self, key, analyte: 'Analyte'):
        self.analytes[str(key)] = analyte
        self._update_mixture_members_term()

    def remove_analyte(self, analyte_id):
        del self.analytes[str(analyte_id)]
        self._update_mixture_members_term()

    def has_analyte(self, analyte_id) -> bool:
        return str(analyte_id) in self.analytes

    def get_member_interpretation(self, member_id) -> 'InterpretationMember':
        return self.member_interpretations[str(member_id)]

    def add_member_interpretation(self, interpretation_member: 'InterpretationMember'):
        self.set_member_interpretation(interpretation_member.id, interpretation_member)

    def set_member_interpretation(self, key, interpretation_member: 'InterpretationMember'):
        self.member_interpretations[str(key)] = interpretation_member

    def remove_member_interpretation(self, member_id):
        del self.member_interpretations[str(member_id)]

    def validate(self) -> bool:
        """
        Perform validation on each component to confirm this object is well formed.

        Returns
        -------
        bool
        """
        analyte_ids = set(self.analytes)
        member_ids = set(self.member_interpretations)
        valid = True
        if not (analyte_ids >= member_ids):
            warnings.warn(
                f"Interpretation has InterpretationMembers {member_ids - analyte_ids} lacking Analytes")
            valid = False
        return valid

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
        a = ''
        if self.attributes:
            a = f', {self.attributes}'
        return f"{self.__class__.__name__}({d}{a})"


class InterpretationMember(IdentifiedAttributeManager):
    __slots__ = ()


class ProteinDescription(AttributeProxy):
    accession = AttributeManagedProperty("MS:1000885|protein accession")
    name = AttributeManagedProperty("MS:1000886|protein name")
    missed_cleavages = AttributeManagedProperty[int]("MS:1003044|number of missed cleavages")
    cleavage_agent = AttributeManagedProperty("MS:1001045|cleavage agent name")
    number_of_enzymatic_termini = AttributeManagedProperty[int]("MS:1003048|number of enzymatic termini")
    flanking_n_terminal_residue = AttributeManagedProperty("MS:1001112|n-terminal flanking residue")
    flanking_c_terminal_residue = AttributeManagedProperty("MS:1001113|c-terminal flanking residue")


class Analyte(IdentifiedAttributeManager):
    __slots__ = ()

    mass = AttributeManagedProperty[float]("MS:1001117|theoretical mass")
    peptide = AttributeManagedProperty[str]("MS:1003169|proforma peptidoform sequence")
    proteins = AttributeGroupFacet[ProteinDescription](ProteinDescription)

    @property
    def charge(self) -> Optional[int]:
        if self.has_attribute(CHARGE_STATE):
            return self.get_attribute(CHARGE_STATE)
        elif self.has_attribute(PROFORMA_ION):
            ion_val = self.get_attribute(PROFORMA_ION)
            val = proforma.ProForma.parse(ion_val)
            return val.charge_state
        else:
            return None

    @charge.setter
    def charge(self, value):
        if value is not None:
            if self.has_attribute(CHARGE_STATE):
                self.replace_attribute(CHARGE_STATE, value)
            else:
                self.add_attribute(CHARGE_STATE, value)
        else:
            self.remove_attribute(CHARGE_STATE)