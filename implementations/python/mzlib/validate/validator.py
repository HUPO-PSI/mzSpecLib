import itertools
import logging
import warnings
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Deque, Dict, Iterator, List, Optional, Sequence, Tuple, Union

from psims.controlled_vocabulary.entity import Entity, ListOfType

from mzlib.attributes import Attribute, Attributed

from mzlib.spectrum import Spectrum
from mzlib.analyte import Analyte, Interpretation
from mzlib.spectrum_library import SpectrumLibrary

from mzlib.ontology import _VocabularyResolverMixin


from mzlib.validate.level import RequirementLevel
from mzlib.validate.semantic_rule import ScopedSemanticRule, load_rule_set
from mzlib.validate.object_rule import ScopedObjectRuleBase, SpectrumPeakAnnotationRule, ValidationWarning
from mzlib.defaults import DEFAULT_UNITS

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def walk_children(term: Entity):
    queue = Deque([term])
    while queue:
        term = queue.popleft()
        yield term
        queue.extend(term.children)


_is_curie = re.compile(r"(([A-Z]+):(.+))(?:\|.+)?")

def is_curie(value: str) -> bool:
    if isinstance(value, str):
        return _is_curie.match(value)
    return False


@dataclass
class ValidationContext:
    attributes_visited: Dict[Tuple[str, str], bool] = field(default_factory=dict)
    rule_states: Dict[str, Any] = field(default_factory=dict)

    def clear_attributes(self):
        self.attributes_visited.clear()

    def record_attribute(self, attribute: Union[Tuple[str, str], Attribute], result: bool):
        if isinstance(attribute, Attribute):
            attribute = (attribute.key, attribute.group_id)
        self.attributes_visited[attribute] = result

    def visited_attribute(self, attribute: Union[Tuple[str, str], Attribute]) -> bool:
        if isinstance(attribute, Attribute):
            attribute = (attribute.key, attribute.group_id)
        return attribute in self.attributes_visited



def _warning_iterator(iterator: Iterator[Spectrum]) -> Iterator[Spectrum]:
    while True:
        try:
            with warnings.catch_warnings(record=True) as w:
                value = next(iterator)
            vw = [a for a in w if issubclass(a.category, ValidationWarning)]
            yield value, vw
        except StopIteration:
            break
        except:
            raise


def _is_of_type(attrib, relation) -> bool:
    if isinstance(relation.value_type.type_definition, type):
        return isinstance(attrib.value, relation.value_type.type_definition)
    else:
        return _try_convert(attrib.value, relation.value_type.type_definition)


def _try_convert(value, converter):
    try:
        converter(value)
        return True
    except (ValueError, TypeError):
        return False


class ValidatorBase(_VocabularyResolverMixin):
    error_log: List
    current_context: ValidationContext

    def reset_context(self):
        self.current_context.clear_attributes()

    def add_warning(self, obj: Attributed, path: str, identifier_path: Tuple, attrib: Any, value: Any, requirement_level: RequirementLevel, message: str):
        raise NotImplementedError()

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary, parsing_warnings: Optional[List[warnings.WarningMessage]] = None):
        raise NotImplementedError()

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        raise NotImplementedError()

    def check_attributes(self, obj: Attributed, path: str, identifer_path: Tuple) -> bool:
        valid: bool = True
        for attrib in obj.attributes:
            if attrib.key == "MS:1003276|other attribute value" or attrib.key == "MS:1003275|other attribute name":
                continue
            if self.current_context.visited_attribute(attrib):
                continue
            try:
                term = self.find_term_for(attrib.key.split("|")[0])
            except KeyError as err:
                logger.warn(f"Could not resolve term for {attrib.key}")
                continue
            value_types = term.get('has_value_type')
            if not value_types:
                value_parsed = isinstance(attrib.value, Attribute)
                if value_parsed or is_curie(attrib.value):
                    if value_parsed:
                        value_key = attrib.value.key.split("|")[0]
                    else:
                        value_key = attrib.value.split("|")[0]
                    if is_curie(value_key):
                        value_term = self.find_term_for(value_key)
                    else:
                        value_term = None
                    if not value_term or not value_term.is_of_type(term):
                        self.add_warning(obj, path, identifer_path, attrib.key, attrib.value, RequirementLevel.must,
                                         f"The value type of {attrib.key} must be a term derived from {attrib.key}, but found {attrib.value}")
                        valid = False
                        continue
                else:
                    self.add_warning(obj, path, identifer_path, attrib.key, attrib.value, RequirementLevel.must,
                                     f"The value type of {attrib.key} must be a term derived from {attrib.key}")
                    valid = False
                    continue
            else:
                for rel in value_types:
                    if isinstance(rel.value_type, ListOfType):
                        hit = False
                        for tp in rel.value_type.type_definition.entity.has_value_type:
                            if isinstance(attrib.value, Sequence) and all(isinstance(v, tp.value_type.type_definition) for v in attrib.value):
                                hit = True
                                break
                        if hit:
                            break
                    elif _is_of_type(attrib, rel):
                        break
                else:
                    self.add_warning(obj, path, identifer_path, attrib.key, attrib.value, RequirementLevel.must,
                                     f"The value type of {attrib.key} must be a value of type {', '.join([rel.value_type.id for rel in value_types])}, but got {type(attrib.value)}")
                    valid = False

            units = term.get('has_units')
            if units:
                if not isinstance(units, list):
                    units = [units]

                if attrib.group_id is not None:
                    try:
                        unit_attrib = obj.get_attribute("UO:0000000|unit", group_identifier=attrib.group_id, raw=True)
                    except KeyError:
                        unit_attrib = None
                        if len(units) == 1:
                            logger.warning(f"{attrib.key}'s unit is missing, defaulting to {units[0]}")
                            continue
                else:
                    unit_attrib = None
                    if len(units) == 1:
                        logger.warning(f"{attrib.key}'s unit is missing, defaulting to {units[0]}")
                        continue
                if unit_attrib:
                    unit_acc, unit_name = unit_attrib.value.split("|", 1)
                    for unit in units:
                        if unit_acc == unit.accession or unit_name == unit.comment:
                            break
                    else:
                        self.add_warning(obj, path, identifer_path, attrib.key, attrib.value, RequirementLevel.must,
                                        f"The attribute {attrib.key} must have a unit {', '.join([rel.accession + '|' + rel.comment for rel in units])}, but got {unit_acc}|{unit_name}")
                        valid = False
                else:
                    if not term.id in DEFAULT_UNITS:
                        self.add_warning(obj, path, identifer_path, attrib.key, attrib.value, RequirementLevel.must,
                                        f"The attribute {attrib.key} must have a unit {', '.join([rel.accession + '|' + rel.comment for rel in units])}, but none were found")
                        valid = False

        return valid

    def validate_library(self, library: SpectrumLibrary, spectrum_iterator: Optional[Iterator[Spectrum]]=None):
        path = "/Library"
        result = self.apply_rules(library, path, (library.identifier, ))
        result &= self.check_attributes(library, path, (library.identifier, ))
        self.reset_context()

        if spectrum_iterator is None:
            spectrum_iterator = library
        for spectrum, warns in _warning_iterator(spectrum_iterator):
            result &= self.validate_spectrum(spectrum, path, library, parsing_warnings=warns)
        return result

    def chain(self, validator: 'ValidatorBase') -> 'ValidatorBase':
        return ValidatorChain([self, validator])

    def __or__(self, other: 'ValidatorBase') -> 'ValidatorBase':
        return self.chain(other)

    def walk_terms_for(self, curie: str) -> Iterator[str]:
        term = self.find_term_for(curie)
        for entity in walk_children(term):
            yield f"{entity.id}|{entity.name}"


@dataclass
class ValidationError:
    path: str
    identifier_path: Tuple
    attribute: Any
    value: Any
    requirement_level: RequirementLevel
    message: str


class Validator(ValidatorBase):
    name: str
    semantic_rules: List[ScopedSemanticRule]
    object_rules: List[ScopedObjectRuleBase]

    def __init__(self, name, semantic_rules: Optional[List[ScopedSemanticRule]] = None, object_rules: Optional[List[ScopedObjectRuleBase]] = None, error_log: Optional[List] = None):
        super().__init__()
        self.name = name
        self.semantic_rules = semantic_rules or []
        self.object_rules = object_rules or []
        self.error_log = error_log or []
        self.current_context = ValidationContext()

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        result = True
        for rule in itertools.chain(self.semantic_rules, self.object_rules):
            if rule.path == path:
                v = rule(obj, path, identifier_path, self)
                level = logging.DEBUG
                if not v and rule.requirement_level > RequirementLevel.may:
                    level = logging.WARN
                    result &= v
                logger.log(level, f"Applied {rule.id} to {path}:{identifier_path} {v}/{result}")
        return result

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary, parsing_warnings: Optional[List[warnings.WarningMessage]] = None):
        path = f"{path}/Spectrum"
        identifier_path = (spectrum.key, )
        result = self.apply_rules(spectrum, path, identifier_path)
        result &= self.check_attributes(spectrum, path, identifier_path)
        self.reset_context()

        if parsing_warnings:
            result = False
            for parsing_warning in parsing_warnings:
                logger.warn(str(parsing_warning.message))

        for _key, analyte in spectrum.analytes.items():
            result &= self.validate_analyte(analyte, path, spectrum, library)

        for _key, interp in spectrum.interpretations.items():
            result &= self.validate_interpretation(interp, path, spectrum, library)
        return result

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        path = f"{path}/Analyte"
        identifier_path = (spectrum.key, analyte.id)
        result = self.apply_rules(analyte, path, identifier_path)
        result &= self.check_attributes(analyte, path, identifier_path)
        self.reset_context()
        return result

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        path = f"{path}/Interpretation"
        identifier_path = (spectrum.key, interpretation.id)
        result = self.apply_rules(interpretation, path, identifier_path)
        result &= self.check_attributes(interpretation, path, identifier_path)
        self.reset_context()
        return result

    def add_warning(self, obj: Attributed, path: str, identifier_path: Tuple, attrib: Any, value: Any, requirement_level: RequirementLevel, message: str):
        if hasattr(obj, "key"):
            key = obj.key
        elif hasattr(obj, "id"):
            key = obj.id
        else:
            key = ''
        warning = f"{attrib.id if hasattr(attrib, 'id') else attrib} failed to validate {path}:{key} ({requirement_level.name.upper()}): {message}"
        logger.warning(warning)
        self.error_log.append(ValidationError(path, identifier_path, attrib, value, requirement_level, warning))

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name!r}, {self.semantic_rules})"


class ValidatorChain(ValidatorBase):
    validators: List[ValidatorBase]

    def __init__(self, validators: List[ValidatorBase], *args, **kwargs):
        self.validators = list(validators)
        super().__init__(*args, **kwargs)

    @property
    def error_log(self):
        log = []
        for validator in self.validators:
            log.extend(validator.error_log)
        return log

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary, parsing_warnings: Optional[List[warnings.WarningMessage]] = None):
        result = True
        for validator in self.validators:
            result &= validator.validate_spectrum(spectrum, path, library, parsing_warnings)
        return result

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        result = True
        for validator in self.validators:
            result &= validator.validate_analyte(analyte, path, spectrum, library)
        return result

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        result = True
        for validator in self.validators:
            result &= validator.validate_interpretation(interpretation, path, spectrum, library)
        return result

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        result = True
        for validator in self.validators:
            result &= validator.apply_rules(obj, path, identifier_path)
        return result

    def check_attributes(self, obj: Attributed, path: str, identifer_path: Tuple) -> bool:
        result = True
        for validator in self.validators:
            result &= validator.check_attributes(obj, path, identifer_path)
        return result

    def reset_context(self):
        for validator in self.validators:
            validator.reset_context()

    def chain(self, validator: ValidatorBase) -> ValidatorBase:
        self.validators.append(validator)
        return self


predicates = {
    "single_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
    "consensus_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
}


object_rules = {
    "peak_annotations": [SpectrumPeakAnnotationRule()]
}


def get_validator_for(name: str) -> Validator:
    rules = load_rule_set(name)
    validator = Validator(name, rules)
    return validator


def get_object_validator_for(name: str) -> Validator:
    rules = object_rules[name]
    validator = Validator(name, object_rules=rules)
    return validator

