import dataclasses
from datetime import datetime
from importlib import resources
import io
import logging

from xml.etree import ElementTree as etree

from typing import Any, List, TYPE_CHECKING, Optional, Tuple, Union

from mzlib.attributes import Attributed
from mzlib.utils import flatten, ensure_iter

from .level import RequirementLevel, CombinationLogic

if TYPE_CHECKING:
    from .validator import ValidatorBase


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class AttributeSemanticPredicate:
    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase") -> bool:
        raise NotImplementedError()


class ValueOfType(AttributeSemanticPredicate):
    type_name: Union[str, List[str]]

    def __init__(self, type_name):
        super().__init__()
        self.type_name = type_name

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase") -> bool:
        if not isinstance(self.type_name, list):
            type_name = [self.type_name]
        else:
            type_name = self.type_name
        result = False
        for tp in type_name:
            if tp == 'int':
                result |= isinstance(value, int) or (isinstance(value, float) and value.is_integer())
            elif tp == 'float':
                result |= isinstance(value, float)
            elif tp == "string":
                result |= isinstance(value, str)
            elif tp == 'datetime':
                if isinstance(value, datetime.datetime):
                    result |= True
                elif isinstance(value, str):
                    logger.warning("Ambiguous datetime found for %s: %s", attribute, value)
                    result |= True
            else:
                raise ValueError(f"Can't validate type {tp}")
        return result


class ValueIsChildOf(AttributeSemanticPredicate):
    parent_accession: str

    def __init__(self, parent_accession):
        super().__init__()
        self.parent_accession = parent_accession

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase"):
        for term in validator_context.walk_terms_for(self.parent_accession):
            if term == value:
                if not term.startswith(self.parent_accession):
                    return True
        return False





@dataclasses.dataclass(frozen=True)
class AttributeSemanticRule:
    accession: str
    name: str
    repeatable: bool
    allow_children: bool
    value: Optional[AttributeSemanticPredicate] = dataclasses.field(default=None)
    condition: Optional[str] = dataclasses.field(default=None)

    @property
    def attribute(self) -> str:
        return f"{self.accession}|{self.name}"


@dataclasses.dataclass
class ScopedSemanticRule:
    id: str
    path: str
    requirement_level: RequirementLevel
    combination_logic: CombinationLogic
    attributes: List[AttributeSemanticRule]

    def find_all_children_of(self, attribute_rule: AttributeSemanticRule, obj: Attributed, validator_context: "ValidatorBase") -> Tuple:
        result = []
        for attrib in validator_context.walk_terms_for(attribute_rule.accession):
            try:
                hits = obj.get_attribute(attrib)
                if hits is not None:
                    result.append((attrib, hits))
            except KeyError:
                continue
        return tuple(result) if result else None

    def validate(self, obj: Attributed, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        values = []

        for attrib in self.attributes:
            if attrib.allow_children:
                values.append(self.find_all_children_of(attrib, obj, validator_context))
            else:
                try:
                    values.append(obj.get_attribute(attrib.attribute))
                except KeyError:
                    values.append(None)

        result = True
        for attrib, value in zip(self.attributes, values):
            # When `get_attribute` matches an attribute multiple times, this will be
            # a list. `find_all_children_of` explicitly returns a tuple to avoid triggering
            # this.
            if isinstance(value, list) and not attrib.repeatable:
                validator_context.add_warning(obj, path, self, identifier_path, value,
                                              self.requirement_level, f"{attrib.attribute} cannot be repeated")
                result = False

        has_value = [v is not None for v in values]
        n_had_values = sum(has_value)

        if self.combination_logic == CombinationLogic.or_:
            if n_had_values < 1:
                validator_context.add_warning(
                    obj, path, identifier_path, self, value, self.requirement_level,
                    f"{self.id} requires one of {', '.join(attr.attribute for attr in self.attributes)}, none found")
                result = False

        elif self.combination_logic == CombinationLogic.and_:
            if n_had_values != len(self.attributes):
                missing = []
                for attr, v in zip(self.attributes, has_value):
                    if not v:
                        missing.append(attr)
                validator_context.add_warning(
                    obj, path, identifier_path, self, value, self.requirement_level,
                    f"{self.id} requires one of {', '.join(attr.attribute for attr in self.attributes)}, "
                    f"{', '.join(m.attribute for m in missing)} were absent")
                result = False

        elif self.combination_logic == CombinationLogic.xor:
            if n_had_values != 1:
                validator_context.add_warning(
                    obj, path, identifier_path, self, value, self.requirement_level,
                    f"{self.id} requires exactly one of {', '.join(attr.attribute for attr in self.attributes)}, "
                    "multiple were found.")
                result = False
        return result

    def __call__(self, obj: Attributed, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        return self.validate(obj, path, identifier_path, validator_context)

    @classmethod
    def from_xml(cls, stream: io.IOBase) -> List['ScopedSemanticRule']:
        tree = etree.parse(stream)
        rules = []
        for node in tree.findall(".//CvMappingRule"):
            term_rules = []
            for term in node.findall("CvTerm"):
                term_rules.append(
                    AttributeSemanticRule(
                        term.attrib['termAccession'],
                        term.attrib['termName'],
                        repeatable=term.attrib['isRepeatable'] == 'true',
                        allow_children=term.attrib['allowChildren'] == 'true'
                    )
                )
            attrib = node.attrib
            rule = ScopedSemanticRule(
                attrib['id'],
                attrib['scopePath'],
                requirement_level=RequirementLevel.from_str(attrib['requirementLevel']),
                combination_logic=CombinationLogic.from_str(attrib['cvTermsCombinationLogic']),
                attributes=term_rules
            )
            rules.append(rule)
        return rules


def load_rule_set(name: str) -> List[ScopedSemanticRule]:
    return ScopedSemanticRule.from_xml(resources.open_binary(__name__.replace(".semantic_rule", '') + '.rules', name + '.xml'))
