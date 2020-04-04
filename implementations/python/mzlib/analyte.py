import textwrap

from mzlib.attributes import AttributeManager


class Analyte(AttributeManager):
    __slots__ = ('id', )

    def __init__(self, id, attributes=None):
        self.id = id
        super(Analyte, self).__init__(attributes)

    def __repr__(self):
        template = f"{self.__class__.__name__}(id={self.id}, "
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[])"
            return template
        template += "[\n%s])" % textwrap.indent(',\n'.join(lines), ' ' * 2)
        return template
