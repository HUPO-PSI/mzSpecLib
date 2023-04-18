import textwrap

from typing import (
    Any, DefaultDict, Iterable,
    Iterator, Optional, Tuple,
    Union, List, Dict,
    Generic, TypeVar, Type
)


T = TypeVar('T')


class Attribute(object):
    __slots__ = ("key", "value", "group_id", "owner_id")
    key: str
    value: Union[str, int, float, 'Attribute', List]
    group_id: Optional[str]
    owner_id: int

    def __init__(self, key, value, group_id=None, owner_id=-1):
        self.key = key
        self.value = value
        self.group_id = group_id
        self.owner_id = owner_id

    def copy(self):
        return self.__class__(self.key, self.value, self.group_id, self.owner_id)

    def __getitem__(self, i):
        if i == 0:
            return self.key
        elif i == 1:
            return self.value
        elif i == 2:
            return self.group_id
        elif i == 3:
            return self.owner_id
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.key
        yield self.value
        if self.group_id:
            yield self.group_id
        yield self.owner_id

    def __len__(self):
        if self.group_id is None:
            return 2
        return 3

    def __eq__(self, other):
        if other is None:
            return False
        return self.key == other.key and self.value == other.value and self.group_id == other.group_id

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return f"{self.__class__.__name__}({self.key}, {self.value}, {self.group_id})"

    def __str__(self):
        base = f"{self.key}={self.value}"
        if self.group_id:
            base = f"[{self.group_id}]" + base
        return base

    def __hash__(self):
        key = (self.key, )
        if self.group_id is not None:
            key = (self.key, self.group_id)
        return hash(key)


class AttributeManager(object):
    """A key-value pair store with optional attribute grouping

    Attributes
    ----------
    attributes: List[List]
        The actual attribute name-value pairs with an optional grouping value
    attribute_dict: Dict
        A mapping from attribute name to indices into :attr:`attributes` and
        the group assignments.
    group_dict: Dict
        A mapping from group identifier to indices into :attr:`attributes`
    group_counter: int
        The number of attribute groups assigned.

    """
    attributes: List[Attribute]
    attribute_dict: Dict
    group_dict: Dict
    group_counter: int


    __slots__ = ('attributes', 'attribute_dict', 'group_dict', 'group_counter')

    def __init__(self, attributes: Iterable = None):
        """

        Parameters
        ----------
        attributes : Iterable[list], optional
            Attribute name-value pairs with an optional grouping value. If omitted,
            the attribute store will be empty.
        """
        self.attributes = []

        # Internal index attributes
        self.attribute_dict = {}
        self.group_dict = {}
        self.group_counter = 1
        if attributes is not None:
            self._from_iterable(attributes)

    def get_next_group_identifier(self) -> str:
        """Retrieve the next un-used attribute group identifier
        and increment the internal counter.

        Returns
        -------
        str
        """

        next_value = self.group_counter
        self.group_counter += 1
        return str(next_value)

    #### Add an attribute to the list and update the lookup tables
    def add_attribute(self, key: str, value, group_identifier: Optional[str] = None):
        """Add an attribute to the list and update the lookup tables

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        items = Attribute(key, value, group_identifier)
        self.attributes.append(items)
        index = len(self.attributes) - 1

        #### If there is already one of these, add it to the lists in attribute_dict
        if key in self.attribute_dict:
            self.attribute_dict[key]["indexes"].append(index)
            if group_identifier is not None:
                self.attribute_dict[key]["groups"].append(group_identifier)

        #### Otherwise, create the entry in attribute_dict
        else:
            if group_identifier is not None:
                self.attribute_dict[key] = {"indexes": [
                    index], "groups": [group_identifier]}
            else:
                self.attribute_dict[key] = {"indexes": [index], "groups": []}

        #### If there is a group_identifier, then update the group_dict
        if group_identifier is not None:
            #### If this group already has one or more entries, add to it
            if group_identifier in self.group_dict:
                self.group_dict[group_identifier].append(index)
            #### Else create and entry for the group_identifier
            else:
                self.group_dict[group_identifier] = [index]

    def add_attribute_group(self, attributes: List[Union[Attribute, Tuple[str, Any]]]):
        group_id = self.get_next_group_identifier()
        for attr in attributes:
            if isinstance(attr, Attribute):
                key = attr.key
                value = attr.value
            else:
                key, value = attr
            self.add_attribute(key, value, group_id)

    def get_attribute(self, key: str, group_identifier: Optional[str] = None,
                      raw: bool = False) -> Union[Any, List[Any], Attribute,
                                                  List[Attribute]]:
        """Get the value or values associated with a given
        attribute key.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                if raw:
                    return [self.attributes[i] for i in indices]
                return [self.attributes[i][1] for i in indices]
            else:
                if raw:
                    return self.attributes[indices[0]]
                return self.attributes[indices[0]][1]
        else:
            groups = indices_and_groups['groups']
            try:
                i = groups.index(group_identifier)
            except ValueError:
                raise KeyError(f"Group {group_identifier} not found") from None
            indices = indices_and_groups['indexes']
            idx = indices[i]
            if raw:
                return self.attributes[idx]
            return self.attributes[idx][1]

    def get_attribute_group(self, group_identifier: str) -> List[Any]:
        result = []
        group_identifier = str(group_identifier)
        for k, indices_and_groups in self.attribute_dict.items():
            for i, g in zip(indices_and_groups['indexes'], indices_and_groups['groups']):
                if g == group_identifier:
                    result.append(self.attributes[i])
        return result

    def replace_attribute(self, key, value, group_identifier=None):
        try:
            indices_and_groups = self.attribute_dict[key]
        except KeyError:
            return self.add_attribute(key, value, group_identifier=group_identifier)
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                raise ValueError("Cannot replace the value of an attribute used multiple times")
            else:
                self.attributes[indices[0]].value = value
        else:
            raise NotImplementedError()

    def get_by_name(self, name: str):
        '''Search for an attribute by human-readable name.

        Parameters
        ----------
        name: str
            The name to search for.

        Returns
        -------
        object:
            The attribute value if found or :const:`None`.
        '''
        matches = []
        for attr in self:
            if attr.key.split("|")[-1] == name:
                matches.append(attr[1])
        n = len(matches)
        if n == 1:
            return matches[0]
        elif n > 1:
            return matches
        return None

    def clear(self):
        """Remove all content from the store.

        """
        self._clear_attributes()

    def remove_attribute(self, key, group_identifier=None):
        """Remove the value or values associated with a given
        attribute key from the store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indexes']
            if len(indices) > 1:
                indices = sorted(indices, reverse=True)
                for i in indices:
                    self.attributes.pop(i)
            else:
                self.attributes.pop(indices[0])
        else:
            groups = indices_and_groups['groups']
            i = groups.index(group_identifier)
            indices = indices_and_groups['indexes']
            idx = indices[i]
            self.attributes.pop(idx)
        attributes = self.attributes
        self.clear()
        self._from_iterable(attributes)

    def _remove_attribute_group(self, group_identifier):
        group_indices = self.group_dict.pop(group_identifier)
        for offset, i in enumerate(sorted(group_indices)):
            self.attributes.pop(i - offset)
        attributes = self.attributes
        self.clear()
        self._from_iterable(attributes)

    def _iter_attribute_groups(self):
        seen = set()
        for group_id, indices in self.group_dict.items():
            yield group_id, [self.attributes[i] for i in indices]
            seen.update(indices)
        acc = []
        for i, attr in enumerate(self.attributes):
            if i in seen:
                continue
            acc.append(attr)
        yield None, acc

    def has_attribute(self, key):
        """Test for the presence of a given attribute

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        return key in self.attribute_dict

    def __contains__(self, key):
        return self.has_attribute(key)

    def __getitem__(self, key):
        return self.get_attribute(key)

    def __setitem__(self, key, value):
        if self.has_attribute(key):
            self.replace_attribute(key, value)
        else:
            self.add_attribute(key, value)

    def keys(self):
        return self.attribute_dict.keys()

    def __eq__(self, other):
        if other is None:
            return False
        self_keys = self.keys()
        other_keys = other.keys()
        if len(self_keys) != len(other_keys):
            return False
        for key in self_keys:
            try:
                self_value = self.get_attribute(key)
                other_value = other.get_attribute(key)
            except KeyError:
                return False
            if self_value != other_value:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return self._count_attributes()

    def __bool__(self):
        return len(self) > 0

    def __iter__(self):
        return self._iter_attributes()

    def _count_attributes(self) -> int:
        return len(self.attributes)

    def _iter_attributes(self) -> Iterator[Attribute]:
        return iter(self.attributes)

    def _clear_attributes(self):
        self.attributes = []

        # Internal index attributes
        self.attribute_dict = {}
        self.group_dict = {}
        self.group_counter = 1

    def _from_iterable(self, attributes):
        mapping = {}
        for attrib in attributes:
            if len(attrib) == 3:
                # this is a grouped attribute
                if attrib[2] in mapping:
                    remap = mapping[attrib[2]]
                else:
                    remap = self.get_next_group_identifier()
                    mapping[attrib[2]] = remap
                self.add_attribute(attrib[0], attrib[1], remap)
            else:
                self.add_attribute(attrib[0], attrib[1])

    def _attributes_from_iterable(self, attributes):
        return self._from_iterable(attributes)

    def copy(self):
        """Make a deep copy of the object
        """
        return self.__class__(self.attributes)

    def __repr__(self):
        if len(self) == 0:
            return f"{self.__class__.__name__}([])"
        lines = list(map(str, self.attributes))
        template = "{}([\n{}])"
        return template.format(
            self.__class__.__name__,
            textwrap.indent(',\n'.join(lines), ' ' * 2))


class IdentifiedAttributeManager(AttributeManager):
    __slots__ = ('id', )

    id: str

    def __init__(self, id, attributes: Iterable = None):
        self.id = str(id)
        super(IdentifiedAttributeManager, self).__init__(attributes)

    def __repr__(self):
        template = f"{self.__class__.__name__}(id={self.id}, "
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[])"
            return template
        template += "[\n%s])" % textwrap.indent(',\n'.join(lines), ' ' * 2)
        return template


class _ReadAttributes(object):
    __slots__ = ()

    attributes: AttributeManager

    def get_attribute(self, key, group_identifier=None, raw: bool = False):
        """Get the value or values associated with a given
        attribute key from the entity's attribute store.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        return self.attributes.get_attribute(key, group_identifier=group_identifier, raw=raw)

    def get_attribute_group(self, group_identifier: str) -> List[Any]:
        return self.attributes.get_attribute_group(group_identifier)

    def has_attribute(self, key) -> bool:
        """Test for the presence of a given attribute in the library
        level store.

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        return self.attributes.has_attribute(key)

    def get_by_name(self, name: str):
        '''Search for an attribute by human-readable name.

        Parameters
        ----------
        name: str
            The name to search for.

        Returns
        -------
        object:
            The attribute value if found or :const:`None`.
        '''
        return self.attributes.get_by_name(name)

    def _iter_attribute_groups(self):
        return self.attributes._iter_attribute_groups()

    def _count_attributes(self) -> int:
        return self.attributes._count_attributes()

    def _iter_attributes(self) -> Iterator[Attribute]:
        return self.attributes._iter_attributes()


class _WriteAttributes(object):
    __slots__ = ()

    attributes: AttributeManager

    def add_attribute(self, key, value, group_identifier=None) -> Union[Any, List[Any]]:
        """Add an attribute to the entity's attributes store.

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        return self.attributes.add_attribute(key, value, group_identifier=group_identifier)

    def replace_attribute(self, key, value, group_identifier=None):
        return self.attributes.replace_attribute(key, value, group_identifier=group_identifier)

    def remove_attribute(self, key, group_identifier=None):
        """Remove the value or values associated with a given
        attribute key from the entity's attribute store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        return self.attributes.remove_attribute(key, group_identifier=group_identifier)

    def _attributes_from_iterable(self, attributes):
        return self.attributes._attributes_from_iterable(attributes)

    def _clear_attributes(self):
        return self.attributes._clear_attributes()


class AttributedEntity(_ReadAttributes, _WriteAttributes):
    '''A base type for entities which contain an :class:`AttributeManager`
    without being completely subsumed by it.

    An :class:`AttributeManager` represents a collection of attributes
    first and foremost, supplying :class:`~.collections.abc.MutableMapping`-like
    interface to them, in addition to methods.
    '''
    __slots__ = ("attributes", )

    attributes: AttributeManager

    def __init__(self, attributes: Iterable=None, **kwargs):
        self.attributes = AttributeManager(attributes)
        super().__init__(**kwargs)


Attributed = Union[AttributeManager, AttributedEntity]


class AttributeManagedProperty(Generic[T]):
    __slots__ = ("attribute", )
    attribute: str

    def __init__(self, attribute: str):
        self.attribute = attribute

    def _get_group_id(self, inst: AttributeManager) -> Optional[str]:
        return getattr(inst, "group_id", None)

    def __get__(self, inst: AttributeManager, cls: Type) -> T:
        if inst is None:
            return self
        if inst.has_attribute(self.attribute):
            return inst.get_attribute(self.attribute, group_identifier=self._get_group_id(inst))
        return None

    def __set__(self, inst: AttributeManager, value: T):
        attrib = self.attribute
        if inst.has_attribute(attrib):
            inst.replace_attribute(attrib, value, group_identifier=self._get_group_id(inst))
        elif inst._count_attributes() > 0:
            attribs = [[attrib, value]] + list(inst._iter_attributes())
            inst._clear_attributes()
            inst._attributes_from_iterable(attribs)
        else:
            inst.add_attribute(attrib, value, group_identifier=self._get_group_id(inst))

    def __delete__(self, inst: AttributeManager, attr):
        inst.remove_attribute(self.attribute, group_identifier=self._get_group_id(inst))


class AttributeListManagedProperty(Generic[T]):
    __slots__ = ("attributes", )
    attributes: List[str]

    def __init__(self, attributes: List[str]):
        self.attributes = attributes

    def __get__(self, inst: AttributeManager, cls: Type) -> T:
        if inst is None:
            return self
        key, val = self._find_key_used(inst)
        if key is None:
            raise KeyError(self.attributes[0])
        return val

    def _find_key_used(self, inst: AttributeManager) -> Optional[Tuple[str, T]]:
        for attr in self.attributes:
            try:
                return attr, inst.get_attribute(attr)
            except KeyError:
                continue
        return None, None

    def __set__(self, inst: AttributeManager, value: T):
        key, _val = self._find_key_used(inst)
        if key is None:
            attrib = self.attributes[0]
        else:
            attrib = key
        if inst.has_attribute(attrib):
            inst.replace_attribute(attrib, value)
        elif inst._count_attributes() > 0:
            attribs = [[attrib, value]] + list(inst._iter_attributes())
            inst._clear_attributes()
            inst._attributes_from_iterable(attribs)
        else:
            inst.add_attribute(attrib, value)

    def __delete__(self, inst: AttributeManager, attr):
        key, val = self._find_key_used(inst)
        if key is not None:
            inst.remove_attribute(key)


class AttributeProxyMeta(type):
    def __new__(cls, typename, bases, namespace):
        props = {}
        for k, v in namespace.items():
            if isinstance(v, (AttributeManagedProperty, AttributeListManagedProperty)):
                props[k] = v
        for base in bases[::-1]:
            if hasattr(base, '_attribute_props'):
                for k, v in base._attribute_props.items():
                    if k not in props:
                        props[k] = v
        namespace['_attribute_props'] = props
        if "__repr__" not in namespace:
            def __repr__(self):
                pairs = []
                for k in self._attribute_props:
                    try:
                        value = getattr(self, k)
                    except KeyError:
                        value = None
                    pairs.append((k, value))
                content = ', '.join(f"{k}={v!r}" for k, v in pairs)
                return f"{self.__class__.__name__}({content})"
            namespace['__repr__'] = __repr__
        return super().__new__(cls, typename, bases, namespace)


class ROAttributeProxy(_ReadAttributes, metaclass=AttributeProxyMeta):
    attributes: Attributed

    def __init__(self, attributes):
        self.attributes = attributes


AttributeProxy = ROAttributeProxy


class AttributeSet(AttributedEntity):
    name: str

    def __init__(self, name: str, attributes: Iterable = None, **kwargs):
        super().__init__(attributes, **kwargs)
        self.name = name

    def member_of(self, target: Attributed) -> bool:
        for attrib in self.attributes:
            if attrib.group_id:
                raise NotImplementedError()
            if not target.has_attribute(attrib.key):
                return False
        return True

    def apply(self, target: Attributed):
        terms_to_remove: List[Tuple[str, Union[Attribute, List[Attribute]]]] = []
        for key in self.attributes.keys():
            terms_to_remove.append((key, target.get_attribute(key, raw=True)))

        group_ids = DefaultDict(int)
        for key, terms in terms_to_remove:
            if isinstance(terms, list):
                for term in terms:
                    if term.group_id:
                        group_ids[term.group_id] += 1
                    target.remove_attribute(key, group_identifier=term.group_id)
            else:
                if term.group_id:
                    group_ids[term.group_id] += 1
                target.remove_attribute(key, group_identifier=term.group_id)

        for group_id in group_ids:
            target._remove_attribute_group(group_id)

        for group_id, attrs in self._iter_attribute_groups():
            if group_id is None:
                for a in attrs:
                    target.add_attribute(a)
            else:
                target.add_attribute_group(attrs)

    def __repr__(self):
        template = f"{self.__class__.__name__}(name={self.name}, "
        lines = list(map(str, self.attributes))
        if not lines:
            template += "[])"
            return template
        template += "[\n%s])" % textwrap.indent(',\n'.join(lines), ' ' * 2)
        return template


class AttributeFacet(Generic[T]):
    facet_type: Type[T]

    def __init__(self, facet_type: Type[T]):
        self.facet_type = facet_type

    def __get__(self, inst: Attributed, klass: Type[Attributed]) -> T:
        if inst is not None:
            return self.facet_type(inst)
        else:
            return self


class AttributeGroupFacet(Generic[T]):
    facet_type: Type[T]

    def __init__(self, facet_type: Type[T]):
        self.facet_type = facet_type

    def wraps(self, attributed: Attributed) -> List[T]:
        attribute_groups = []
        facet_keys = {a.attribute for a in self.facet_type._attribute_props.values()}
        for group_id in attributed.group_dict:
            group = attributed.get_attribute_group(group_id)
            group_keys = {a.key for a in group}
            if group_keys & facet_keys:
                attribute_groups.append(group_id)

        facets = []
        for group_id in attribute_groups:
            facet = self.facet_type(attributed)
            facet.group_id = group_id
            facets.append(facet)
        return facets

    def __get__(self, inst: Attributed, klass) -> Union[List[T], 'AttributeGroupFacet']:
        if inst is not None:
            return self.wraps(inst)
        return self

