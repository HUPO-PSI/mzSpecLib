import logging
from typing import Any, Deque, Dict, Iterator, List, Optional, Tuple

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
    def add_warning(self, obj: Attributed, path: str, identifier_path: Tuple, attrib: Any, value: Any, requirement_level: RequirementLevel, message: str):
        raise NotImplementedError()

    def validate_spectrum(self, spectrum: Spectrum, path: str, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_analyte(self, analyte: Analyte, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_interpretation(self, interpretation: Interpretation, path: str, spectrum: Spectrum, library: SpectrumLibrary):
        raise NotImplementedError()

    def validate_library(self, library: SpectrumLibrary):
        for spectrum in library:
            self.validate_spectrum(spectrum, "/Library")

    def walk_terms_for(self, curie: str) -> Iterator[str]:
        term = self.find_term_for(curie)
        for entity in walk_children(term):
            yield f"{entity.id}|{entity.name}"


class Validator(ValidatorBase):
    name: str
    semantic_rules: List[ScopedSemanticRule]
    object_rules: List
    error_log: List

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
                logger.info(f"Applied {rule.id} to {path}:{identifier_path} {v}/{result}")
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
        warning = f"{attrib} failed to validate {path}:{key} ({requirement_level.name.upper()}): {message}"
        logger.warn(warning)
        self.error_log.append(warning)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name!r}, {self.semantic_rules})"


def get_validator_for(name: str) -> Validator:
    rules = load_rule_set(name)
    return Validator(name, rules)