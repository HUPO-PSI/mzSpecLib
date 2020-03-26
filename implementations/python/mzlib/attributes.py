
class AttributeManager(object):

    def __init__(self, attributes=None):
        self.attributes = []

        # Internal index attributes
        self.attribute_dict = {}
        self.group_dict = {}
        self.group_counter = 1
        if attributes is not None:
            self._from_iterable(attributes)

    #### Get the next group identifier
    def get_next_group_identifier(self):
        next = self.group_counter
        self.group_counter += 1
        return(str(next))

    #### Add an attribute to the list and update the lookup tables
    def add_attribute(self, key, value, group_identifier=None):
        items = [key, value]
        if group_identifier is not None:
            items.append(group_identifier)
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

    def get_attribute(self, key, group_identifier=None):
        indices_and_groups = self.attribute_dict[key]
        if group_identifier is None:
            indices = indices_and_groups['indices']
            if len(indices) > 1:
                return [self.attributes[i] for i in indices]
            else:
                return self.attributes[indices[0]]
        else:
            groups = indices_and_groups['groups']
            i = groups.index(group_identifier)
            indices = indices['indices']
            idx = indices[i]
            return self.attributes[idx]

    def _from_iterable(self, attributes):
        for attrib in attributes:
            if len(attrib) == 3:
                # this is a grouped attribute
                self.add_attribute(attrib[0], attrib[1], attrib[2])
                spec.group_counter = int(attrib[2])
            else:
                self.add_attribute(attrib[0], attrib[1])

    def copy(self):
        return self.__class__(self.attributes)