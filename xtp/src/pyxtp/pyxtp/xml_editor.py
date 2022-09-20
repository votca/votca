"""Module to edit the options in the XML file."""

import xml.etree.ElementTree as ET
from typing import Any, Dict
from types import SimpleNamespace
import xmltodict 
import os

def edit_xml(root: ET.Element, sections: Dict[str, Any], path: str = ".dftgwbse"):
    """Edit a XML object using sections."""
    sec = root.find(path)
    if sec is None:
        raise RuntimeError(f"Unknown Section: {path}")
    else:
        for key, val in sections.items():
            update_node(sec, key, val)


def update_node(root: ET.Element, key: str, value: Any):
    """Update nodes recursively."""
    sec = root.find(key)

    # insert new empty node
    if sec is None:
        #elem = ET.Element(key)
        #root.insert(0, elem)
        #update_node(root, key, value)
        raise RuntimeError(f"Unknown option: {key}")

    else:
        for node in root.findall(key):
            if not isinstance(value, dict):
                node.text = str(value)
            else:
                for k, v in value.items():
                    update_node(node, k, v)



def create_xml_tree(root: ET.Element, dict_tree: Dict[str, Any] ):
    # Node : recursively create tree nodes
    if type(dict_tree) == dict:
        for k, v in dict_tree.items():
            create_xml_tree(ET.SubElement(root, k), v)
        return root
    # Leaf : just set the value of the current node
    else:
        root.text = str(dict_tree)
        
class NestedNamespace(SimpleNamespace):
    """Extend SimpleNamespace to nested structure

    example:
    .. code-block:: python
       d = NestedNamespace({'a': 1, 'b':{'x':2, 'y':3}})
       d.a    # 1
       d.b.x  # 2

    In case of elements starting with a @ they are replaced
    by underscore, e.g. :

    .. code-block:: python
       d = NestedNamespace({'a': 1, 'b':{'@x':2, 'y':3}})
       d.a    # 1
       d.b._x  # 2
    """

    def __init__(self, dictionary, **kwargs):

        super().__init__(**kwargs)

        for key, value in dictionary.items():
            if isinstance(value, dict):
                self.__setattr__(key, NestedNamespace(value))
            else:
                if key.startswith('@'):
                    key = '_'+key[1:]
                self.__setattr__(key, value)

def xml2namespace(xml_filename: str, xml_attribs=True) -> NestedNamespace:
    """Load an xml file into an (nested) namesapce using xml2dict

    Args:
        xml_filename (str): name of the XML file
        xml_attribs (bool): load the xml attirubtes

    Returns:
        NestedNamespace: namespace
    """
    with open(xml_filename) as fd:
        dict_data = xmltodict.parse(fd.read(), xml_attribs=xml_attribs)
    return NestedNamespace(dict_data)   

def namespace2dict(namespace_data: NestedNamespace, 
                   path: str='') -> dict:
    """Transform a nested namespace into a dictionary.

    Args:
        namespace_data (NestedNamespace): data in namespace form
        path (str, optional): basepath in the structure. Defaults to ''.
        remove_optional(bool, optional): remove the attr starting with '_' Default to True
    
    Returns:
        dict: dictionary maps the namespace
    """
    output = dict()
    for key, value in namespace_data.__dict__.items():
        if isinstance(value, NestedNamespace):
            output.update(namespace2dict(value, path=os.path.join(path,key)))
        else:
            if not key.startswith('_'):
                output[os.path.join(path,key)] = value
    return output
