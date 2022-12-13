
import xml.etree.ElementTree as ET
from typing import Any, Dict
from types import SimpleNamespace
import xmltodict 
import os
from lxml import etree
from lxml.etree import Element, ElementTree

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
                
    @classmethod
    def __fromxml__(cls, xml_str: str, xml_attribs=True):
        """Creates a class instance from an xml file stored in a string

        Args:
            xml_str (str): content of an xml file represented as string

        Returns:
            NestedNamespace: namesapce mapping the xml file
            
        examle:
        .. code-block:: python
            from lxml import etree
            tree = etree.parse('test.xml')
            xml_str = etree.tostring(tree)
            data = NestedNamespace._fromxml(xml_str)
        """
        dict_data = xmltodict.parse(xml_str, xml_attribs=xml_attribs)
        return cls(dict_data)   
    
    def __todict__(self) -> dict:
        """Transform the nested namespace into a dictionary.
        
        Returns:
            dict: dictionary maps the namespace
        """
        
        def _recursive_namespace2dict(namespace_data: NestedNamespace, 
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
                    output.update(_recursive_namespace2dict(value, path=os.path.join(path,key)))
                else:
                    if not key.startswith('_'):
                        output[os.path.join(path,key)] = value
            return output

        return _recursive_namespace2dict(self)
    
    def _haschildren(self) -> bool:
        """Returns true if self contains other nestednamespaces

        Returns:
            bool: _description_
        """
        for k,v in self.__dict__.items():
            if isinstance(v, NestedNamespace):
                return True
        return False