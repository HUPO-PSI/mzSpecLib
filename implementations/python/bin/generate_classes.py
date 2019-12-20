#!/usr/bin/env python3

import os
import yaml


PROJECT_PATH = 'mzlib'


def main():
    base_path = os.path.dirname(os.path.abspath(__file__))
    classes_yaml_file = base_path + "/../etc/classes.yaml"
    print(f"INFO: Reading from {classes_yaml_file}")

    with open(classes_yaml_file, 'r') as stream:
        definitions = yaml.safe_load(stream)

    try:
        classes = definitions["components"]["schemas"]
    except:
        print(f"ERROR: Did not find [components][schemas] in {classes_yaml_file}")
        return()

    for class_name in classes:
        print(f"INFO: Generating code for class {class_name}")
        class_description = classes[class_name]["description"] if "description" in classes[class_name] else ""
        with open(f"{base_path}/../{PROJECT_PATH}/{class_name}.template.py", 'w') as outfile:
            outfile.write(f'''#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


''')

            properties_def_buffer = ""
            properties_series_buffer = ""
            properties_doc_buffer = ""
            class_doc_buffer = f'''    """
    {class_name} - {class_description}

    Attributes
    ----------'''
            constructor_buffer = ""
            methods_def_buffer = ""
            methods_doc_buffer = ""
            attributes_doc_buffer = ""
            init_setters_buffer = ""

            if "properties" in classes[class_name]:
                for property_name in classes[class_name]["properties"]:
                    property = classes[class_name]["properties"][property_name]
                    description = property["description"] if "description" in property else ""
                    properties_def_buffer += f'''    #### Define getter/setter for attribute {property_name}
    @property
    def {property_name}(self):
        return(self._{property_name})
    @{property_name}.setter
    def {property_name}(self, {property_name}):
        self._{property_name} = {property_name}

'''
                    default = property["default"] if "default" in property else "None"
                    if default != "": default = f"={default}"
                    type = property["type"] if "type" in property else "string"
                    properties_series_buffer += f''', {property_name}{default}'''
                    properties_doc_buffer = f'''        {property_name} : {type}
            {description}
'''
                    attributes_doc_buffer = f'''    {property_name} : {type}
        {description}
'''
                    init_setters_buffer += f'''        self.{property_name} = {property_name}
'''

            constructor_buffer += f'''
    def __init__(self{properties_series_buffer}):
        """
        __init__ - {class_name} constructor

        Parameters
        ----------
{properties_doc_buffer}
        """

{init_setters_buffer}
'''

            if "methods" in classes[class_name]:
                for method_name in classes[class_name]["methods"]:
                    method = classes[class_name]["methods"][method_name]
                    if method is None: continue
                    method_description = method["description"] if "description" in method else ""

                    parameters_doc_buffer = ""
                    parameters_series_buffer = ""
                    if "parameters" in method:
                        for parameter_name in method["parameters"]:
                            parameter = method["parameters"][parameter_name]
                            if parameter is None: continue
                            description = parameter["description"] if "description" in parameter else ""
                            type = parameter["type"] if "type" in parameter else "string"
                            default = parameter["default"] if "default" in parameter else ""
                            if default != "": default = f"={default}"
                            parameters_doc_buffer += f'''        {parameter_name} : {type}
            {description}
'''
                            parameters_series_buffer += f''', {parameter_name}{default}'''

                    methods_def_buffer += f'''    def {method_name}(self{parameters_series_buffer}):
        """
        {method_name} - {method_description}

        Extended description of function.

        Parameters
        ----------
{parameters_doc_buffer}
        Returns
        -------
        int
            Description of return value
        """

        #### Begin functionality here

        return()


'''
                    methods_doc_buffer += f'''    {method_name} - {method_description}
'''                    

            outfile.write(f'''class {class_name}:
{class_doc_buffer}
{attributes_doc_buffer}
    Methods
    -------
{methods_doc_buffer}
    """

{constructor_buffer}
{properties_def_buffer}

{methods_def_buffer}

''')

            outfile.write(f'''

#### Example using this class
def example():

    #### Create a new RTXFeedback object
    obj = {class_name}()
    #obj.dosomething()
    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example()
    return()


if __name__ == "__main__": main()

''')

if __name__ == "__main__": main()
