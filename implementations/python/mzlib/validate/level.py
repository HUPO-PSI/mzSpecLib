import enum

class RequirementLevel(enum.Enum):
    may = enum.auto()
    should = enum.auto()
    must = enum.auto()
    disallow = enum.auto()

    @classmethod
    def from_str(cls, text: str) -> 'RequirementLevel':
        text = str(text).lower()
        return cls[text]


class CombinationLogic(enum.Enum):
    or_ = enum.auto()
    and_ = enum.auto()
    xor = enum.auto()

    @classmethod
    def from_str(cls, text: str) -> 'CombinationLogic':
        text = str(text).lower()
        if text in ("or", "and"):
            text += '_'
        return cls[text]
