import enum

class RequirementLevel(enum.IntEnum):
    may = 0
    should = 1
    must = 2
    disallow = 3

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

    def to_str(self) -> str:
        if self == self.or_:
            return "OR"
        elif self == self.and_:
            return "AND"
        elif self == self.xor:
            return "XOR"
        else:
            raise ValueError(self)
