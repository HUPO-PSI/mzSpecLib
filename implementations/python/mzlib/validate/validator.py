from dataclasses import dataclass
import itertools
import logging
from typing import Any, Callable, Deque, Dict, Iterator, List, Optional, Tuple

from psims.controlled_vocabulary import Entity

from mzlib.attributes import Attributed
from mzlib.spectrum import Spectrum
from mzlib.analyte import Analyte, Interpretation
from mzlib.spectrum_library import SpectrumLibrary

from mzlib.backends.base import VocabularyResolverMixin


from mzlib.validate.level import RequirementLevel
from mzlib.validate.semantic_rule import ScopedSemanticRule, load_rule_set

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def walk_children(term: Entity):
    queue = Deque([term])
    while queue:
        term = queue.popleft()
        yield term
        queue.extend(term.children)


class ValidatorBase(VocabularyResolverMixin):
    error_log: List

    def add_warning(self, obj: Attributed, path: str, identifier_path: Tuple, attrib: Any, value: Any, requirement_level: RequirementLevel, message: str):
        raise NotImplementedError()

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        raise NotImplementedError()

    def validate_library(self, library: SpectrumLibrary, spectrum_iterator: Optional[Iterator[Spectrum]]=None):
        path = "/Library"
        self.apply_rules(library, path, (library.identifier, ))
        result = True
        if spectrum_iterator is None:
            spectrum_iterator = library
        for spectrum in spectrum_iterator:
            result &= self.validate_spectrum(spectrum, path, library)
        return result

    def chain(self, validator: 'ValidatorBase') -> 'ValidatorBase':
        return ValidatorChain([self, validator])

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
    object_rules: List

    def __init__(self, name, semantic_rules: Optional[List[ScopedSemanticRule]]=None, object_rules: Optional[List]=None, error_log: Optional[List]=None):
        super().__init__()
        self.name = name
        self.semantic_rules = semantic_rules or []
        self.object_rules = object_rules or []
        self.error_log = error_log or []

    def apply_rules(self, obj: Attributed, path: str, identifier_path: Tuple) -> bool:
        result = True
        for rule in self.semantic_rules:
            if rule.path == path:
                v = rule(obj, path, identifier_path, self)
                result &= v
                level = logging.DEBUG
                if not v and rule.requirement_level > RequirementLevel.may:
                    level = logging.WARN
                logger.log(level, f"Applied {rule.id} to {path}:{identifier_path} {v}/{result}")
        return result

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary):
        path = f"{path}/Spectrum"
        identifier_path = (spectrum.key, )
        result = self.apply_rules(spectrum, path, identifier_path)

        for _key, analyte in spectrum.analytes.items():
            result &= self.validate_analyte(analyte, path, spectrum, library)

        for _key, interp in spectrum.interpretations.items():
            result &= self.validate_interpretation(interp, path, spectrum, library)
        return result

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        path = f"{path}/Analyte"
        identifier_path = (spectrum.key, analyte.id)
        return self.apply_rules(analyte, path, identifier_path)

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        path = f"{path}/Interpretation"
        identifier_path = (spectrum.key, interpretation.id)
        return self.apply_rules(interpretation, path, identifier_path)

    def add_warning(self, obj: Attributed, path: str, identifier_path: Tuple, attrib: Any, value: Any, requirement_level: RequirementLevel, message: str):
        if hasattr(obj, "key"):
            key = obj.key
        elif hasattr(obj, "id"):
            key = obj.id
        else:
            key = ''
        warning = f"{attrib.id} failed to validate {path}:{key} ({requirement_level.name.upper()}): {message}"
        logger.warning(warning)
        self.error_log.append(ValidationError(path, identifier_path, attrib, value, requirement_level, warning))

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name!r}, {self.semantic_rules})"


class PredicateValidator(Validator):
    spectrum_predicate: Callable[[Spectrum], bool]

    def __init__(self, name, predicate: Callable[[Spectrum], bool], semantic_rules: Optional[List[ScopedSemanticRule]] = None, object_rules: Optional[List] = None, error_log: Optional[List] = None):
        self.spectrum_predicate = predicate
        super().__init__(name, semantic_rules, object_rules, error_log)

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary):
        if self.spectrum_predicate(spectrum):
            return super().validate_spectrum(spectrum, path, library)
        return True


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

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary):
        result = True
        for validator in self.validators:
            result &= validator.validate_spectrum(spectrum, path, library)
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

    def chain(self, validator: ValidatorBase) -> ValidatorBase:
        self.validators.append(validator)
        return self


predicates = {
    "single_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
    "consensus_spectrum": lambda spec: spec.get_attribute("MS:1003065|spectrum aggregation type") == "MS:1003066|singleton spectrum",
}


def get_validator_for(name: str) -> Validator:
    rules = load_rule_set(name)
    return Validator(name, rules)
