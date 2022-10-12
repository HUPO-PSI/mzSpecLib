import dataclasses
import io
import json
import logging
import re

from datetime import datetime
from importlib import resources

from xml.etree import ElementTree as etree

from typing import Any, ClassVar, Dict, List, TYPE_CHECKING, Mapping, Optional, Sequence, Set, Tuple, Union

from mzlib.attributes import Attributed
from mzlib.utils import flatten, ensure_iter
from mzlib.backends.base import VocabularyResolverMixin

from .level import RequirementLevel, CombinationLogic

if TYPE_CHECKING:
    from .validator import ValidatorBase


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class AttributeSemanticPredicate:
    name: ClassVar[str]

    _registry = {}

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase") -> bool:
        raise NotImplementedError()

    def to_dict(self) -> Dict[str, Any]:
        raise NotImplementedError()

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        if isinstance(state, Mapping):
            name = state['name']
        elif isinstance(state, str):
            name = state
        else:
            raise TypeError(f"Cannot convert {state} to {cls.__name__}")
        rule_tp = cls._registry[name]
        return rule_tp.from_dict(state)

    def __init_subclass__(cls, *args, **kwargs):
        cls._registry[cls.name] = cls
        super().__init_subclass__(**kwargs)


class ValueOfType(AttributeSemanticPredicate):
    type_name: Union[str, List[str]]

    name = "value_of_type"

    def __init__(self, type_name):
        super().__init__()
        self.type_name = type_name

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase") -> bool:
        if isinstance(value, list) and attribute.repeatable:
            if all(self.validate(attribute, v, validator_context) for v in value):
                return True
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

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "type_name": self.type_name
        }

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        return cls(state['type_name'])


class ValueIsChildOf(AttributeSemanticPredicate):
    accession: str

    name = "value_is_child_of"

    def __init__(self, accession):
        super().__init__()
        self.accession = accession

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase"):
        if isinstance(value, list) and attribute.repeatable:
            return all(self.validate(attribute, v, validator_context) for v in value)
        for term in validator_context.walk_terms_for(self.accession):
            if term == value:
                if not term.startswith(self.accession):
                    return True
        return False

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "accession": self.accession
        }

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        return cls(state['accession'])


class ValueIsUnique(AttributeSemanticPredicate):
    seen: set

    name = "value_is_unique"

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.seen = set()

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase"):
        if isinstance(value, list) and attribute.repeatable:
            return all(self.validate(attribute, v, validator_context) for v in value)
        if value in self.seen:
            return False
        self.seen.add(value)
        return True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
        }

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        return cls()


class ValueMatches(AttributeSemanticPredicate):
    accession: str

    name = "value_matches"

    def __init__(self, accession):
        super().__init__()
        self.accession = accession

    def validate(self, attribute: 'AttributeSemanticRule', value: str, validator_context: "ValidatorBase"):
        if isinstance(value, list) and attribute.repeatable:
            return all(self.validate(attribute, v, validator_context) for v in value)
        entity = validator_context.find_term_for(self.accession)
        key = f"{entity.id}|{entity.name}"
        return key == value

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "accession": self.accession
        }

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        return cls(state['accession'])


class ValueMatchesPattern(AttributeSemanticPredicate):
    pattern: re.Pattern

    name = "value_matches_pattern"

    def __init__(self, pattern: str):
        super().__init__()
        self.pattern = re.compile(pattern)

    def to_dict(self) -> Dict[str, Any]:
        state = {}
        state['name'] = self.name
        state['pattern'] = self.pattern.pattern
        return state

    @classmethod
    def from_dict(cls, state: Dict[str, Any]) -> 'AttributeSemanticPredicate':
        return cls(state['pattern'])


@dataclasses.dataclass(frozen=True)
class AttributeSemanticRule:
    accession: str
    name: str
    repeatable: bool
    allow_children: bool
    value: Optional[AttributeSemanticPredicate] = dataclasses.field(default=None)
    condition: Optional['AttributeSemanticRule'] = dataclasses.field(default=None)

    @property
    def attribute(self) -> str:
        return f"{self.accession}|{self.name}"

    def to_dict(self) -> Dict[str, Any]:
        state = {
            "accession": self.accession,
            "name": self.name,
            "repeatable": self.repeatable,
            "allow_children": self.allow_children,
        }
        if self.value is not None:
            state['value'] = self.value.to_dict()
        if self.condition is not None:
            state['condition'] = self.condition.to_dict()
        return state

    @classmethod
    def from_dict(cls, state, cv_provider: 'VocabularyResolverMixin') -> 'AttributeSemanticRule':
        if isinstance(state, str):
            state = {
                "accession": state,
            }

        if "name" not in state:
            term = cv_provider.find_term_for(state)
            state['name'] = term.name

        repeatable = state.get("repeatable", False)
        allow_children = state.get("allow_children", False)

        value_rule = state.get("value")
        if value_rule:
            value_rule = AttributeSemanticPredicate.from_dict(value_rule)

        condition = state.get("condition")
        if condition:
            condition = cls.from_dict(condition, cv_provider)

        attr_rule = cls(
            state["accession"],
            state["name"],
            repeatable=repeatable,
            allow_children=allow_children,
            value=value_rule,
            condition=condition,
        )
        return attr_rule


@dataclasses.dataclass
class ScopedSemanticRule:
    id: str
    path: str
    attributes: List[AttributeSemanticRule]
    requirement_level: RequirementLevel
    combination_logic: CombinationLogic
    condition: Optional[AttributeSemanticRule] = dataclasses.field(default=None)

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

    def check_rule(self, obj: Attributed, attrib: AttributeSemanticRule, validator_context: "ValidatorBase") -> bool:
        result = True
        if attrib.allow_children:
            value = self.find_all_children_of(attrib, obj, validator_context)
        else:
            try:
                value = obj.get_attribute(attrib.attribute)
            except KeyError:
                value = None

        if isinstance(value, list) and not attrib.repeatable:
            result = False
        if attrib.value:
            if not attrib.value.validate(attrib, value, validator_context):
                result = False
        return result

    def validate(self, obj: Attributed, path: str, identifier_path: Tuple, validator_context: "ValidatorBase") -> bool:
        values = []

        if self.condition:
            if not self.check_rule(obj, self.condition, validator_context):
                return True

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
                validator_context.add_warning(
                    obj, path, self, identifier_path, value,
                    self.requirement_level,
                    f"{attrib.attribute} cannot be repeated")

                result = False
            if attrib.value:
                if not attrib.value.validate(attrib, value, validator_context):
                    validator_context.add_warning(
                        obj, path, identifier_path, self, value, self.requirement_level,
                        f"{self.id} required the value of {attrib.accession}|{attrib.name} conform to {attrib.value.name}")
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
                id=attrib['id'],
                path=attrib['scopePath'],
                requirement_level=RequirementLevel.from_str(attrib['requirementLevel']),
                combination_logic=CombinationLogic.from_str(attrib['cvTermsCombinationLogic']),
                attributes=term_rules
            )
            rules.append(rule)
        return rules

    @classmethod
    def from_dict(cls, data: Dict[str, Any], cv_provider: 'VocabularyResolverMixin') -> List['ScopedSemanticRule']:
        rules = []
        for rule_spec in data['rules']:
            rule_id = rule_spec['id']
            level = RequirementLevel.from_str(rule_spec['level'].upper())
            path = rule_spec['path']
            combinator = CombinationLogic.from_str(
                rule_spec.get("combination_logic", "OR").upper())

            attribute_rules = []
            for attrib in rule_spec['attr']:
                attr_rule = AttributeSemanticRule.from_dict(attrib, cv_provider)
                attribute_rules.append(attr_rule)

            condition = rule_spec.get("condition")
            if condition:
                condition = AttributeSemanticRule.from_dict(condition, cv_provider)

            rules.append(
                ScopedSemanticRule(
                    id=rule_id,
                    path=path,
                    attributes=attribute_rules,
                    requirement_level=level,
                    combination_logic=combinator,
                    condition=condition
                )
            )
        return rules

    def to_dict(self) -> Dict[str, Any]:
        state = {
            "id": self.id,
            "path": self.path,
            "level": self.requirement_level.name.upper(),
            "combination_logic": self.combination_logic.to_str(),
            "attr": [
                a.to_dict() for a in self.attributes
            ]
        }
        if self.condition:
            state['condition'] = self.condition
        return state


class RuleSet(Sequence[ScopedSemanticRule]):
    name: str
    rules: List[ScopedSemanticRule]

    def __init__(self, name: str, rules: List[ScopedSemanticRule]):
        super().__init__()
        self.name = name
        self.rules = rules

    def __iter__(self):
        return iter(self.rules)

    def __len__(self):
        return len(self.rules)

    def __getitem__(self, i):
        return self.rules[i]

    @classmethod
    def from_dict(cls, data: Dict[str, any]) -> 'RuleSet':
        name = data['name']
        rules = ScopedSemanticRule.from_dict(data)
        return cls(name, rules)

    def to_dict(self) -> Dict[str, Any]:
        state = {
            "name": self.name,
            "rule": [
                rule.to_dict() for rule in self.rules
            ]
        }
        return state


def load_rule_set(name: str) -> List[ScopedSemanticRule]:
    return RuleSet(
        name,
        # ScopedSemanticRule.from_xml(
        #     resources.open_binary(
        #         __name__.replace(".semantic_rule", '') + '.rules',
        #         name + '.xml'
        #     )
        # )
        ScopedSemanticRule.from_dict(
            json.load(resources.open_text(
                __name__.replace(".semantic_rule", '') + '.rules',
                name + '.json'
            )
            ),
            VocabularyResolverMixin()
        )
    )
