"""Module defining the Option class."""
from genericpath import isfile
from multiprocessing.sharedctypes import Value
from typing import Any, Dict
import os
import xml.etree.ElementTree as ET
from lxml import etree
from lxml.etree import Element, ElementTree
from pyxtp.xml_editor import NestedNamespace, xmlstr2namespace, namespace2dict

class Options(dict):
    """Extend the base class dictionary with a '.' notation.

    example:
    .. code-block:: python
       d = Options({'a': 1})
       d['a'] # 1
       d.a    # 1
       d.b = 3
       d["b"] == 3  # True
    """

    def __init__(self, *args, **kwargs):
        """Create a recursive Options object."""
        super().__init__(*args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                self[k] = Options(v)

    def __getattr__(self, key: str):
        """Allow `obj.key` notation."""
        if self.get(key, None) is None:
            self[key] = Options()
        return self.get(key)

    def __setattr__(self, key: str, value: Any):
        """Allow `obj.key = new_value` notation."""
        self.__setitem__(key, value)

    def __deepcopy__(self, _):
        """Deepcopy of the Settings object."""
        return Options(self.copy())

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a normal dictionary."""
        def converter(var):
            return var.to_dict() if (isinstance(var, Options) and len(var) > 0) else var

        return {k: converter(v) for k, v in self.items()}


class XTPOptions(NestedNamespace):

    def __init__(self):
        """Creates a nested namespace data structure from the default input xml files 
        stored in the VOTCASHARE

        Raises:
            ValueError: if VOTCASHARE is not declared as an environement variable
        """

        self.votcashare = os.environ.get('VOTCASHARE')
        if self.votcashare is None:
            raise ValueError('Please set VOTCASHARE environement variable')

        # main filename
        self.dftgwbse_options_file = f'{self.votcashare}/xtp/xml/dftgwbse.xml'      
        
        # parse xmlfile and replace links to the embedded xml files
        self.xml_tree = self.process_xmlfile(self.dftgwbse_options_file)
        
        # create a namespace from the xml tree
        self.data = xmlstr2namespace(etree.tostring(self.xml_tree)).options.dftgwbse
    
    def process_xmlfile(self, xml_filename: str) -> None:
        """Parse the initial xml files and add links

        Args:
            xml_filename (str): filename of the main xmlfile
        """
        
        def _recursive_insert_links(el: Element):
            """recursively add links files to the main tree

            Args:
                el (Element): lxml element

            Raises:
                ValueError: If one of the link is not found
            """
            
            for child in el.getchildren():
                _recursive_insert_links(child)
                
            if 'link' in el.attrib:
                link_name = os.path.join(f'{self.votcashare}/xtp/xml/', el.attrib['link'])
                if not os.path.isfile(link_name):
                    raise ValueError('File %s not found' %link_name)
                xml_link = etree.parse(link_name)
                el.getparent().replace(el, xml_link.getroot())
        
        # parse the file
        tree = etree.parse(xml_filename)
        
        # recursively insert the links
        _recursive_insert_links(tree.getroot())
        
        return tree
    
    def write_xml(self, filename: str) -> None:
        """update and clean the xml and then write the input file

        Args:
            filename (str): filename in the 
        """
        
        # create a dict of the user providedoptions
        dftgwbse_options = namespace2dict(self.data)
        
        # update the options from the user provided input
        self.xml_tree = self._update(self.xml_tree, dftgwbse_options)

        # remove the empty elements
        self.xml_tree = self._clean(self.xml_tree)
        
        # write the file
        self.xml_tree.write(filename)
    

    @staticmethod
    def _update(xml: ElementTree, dict_options: dict) -> ElementTree:
        """Read and update default option files with the user define options

        Args:
            xml_tree (ElementTree): tree of the xml file
            dict_options (dict): dictionary of user defined options

        Returns:
            ET.ElementTree: xml structure containing the update options.
        """

        # replace text value of the dict elements
        for key, value in dict_options.items():
            try:
                xml.find('.//'+key).text = value
            except:
                print(key, ' not found ')
        return xml 
    
    @staticmethod
    def _clean(tree: ET.ElementTree, remove_attributes=[]) -> ElementTree:
        """Recursively remove all empty nodes and elements containing only empty nodes

        Args:
            xml_tree (ElementTree): tree of the xml file
            remove_attributes (list): list of xml attributes to remove
        """
        
        def is_empty(elem: Element) -> bool:
            """returns true if elem does not have text

            Args:
                elem (Element): xml element

            Returns:
                bool true if e is empty
            """
            return (elem.text is None) or (elem.text.strip() == '')
        
        def recursively_remove_empty(elem: Element):
            """reccursive function that remove empty elements

            Args:
                e (Element): xml element

            Returns:
                bool: True if we need to remove the element False otherwise
            """
                   
            # recursive call that gets a list of :
            # [] if e doens't have children
            # [True, False, False, True ] where each bool indicate if e children are empty or not
            child_bool = list(recursively_remove_empty(c) for c in elem.iterchildren())   
            nchildren = len(child_bool)
            
            # there is an issue with the last element
            # so we need to catch the exception
            try:
                
                # if e doesn't have kids, i.e. it's a leaf
                # remove e if it doesn't have text and returns True
                # keep e if it has text and return False
                if nchildren == 0:
                    if is_empty(elem):
                        elem.getparent().remove(elem)
                        return True
                    else:
                        return False
                # remove e if all its children have no text
                # and returns True
                elif all(child_bool):
                    elem.getparent().remove(elem)
                    return True
                
            # excpetion for root ...
            except:        
                return False
            
        
        # remove attributes if any are passed
        etree.strip_attributes(tree, remove_attributes)
        
        # remove the empties
        recursively_remove_empty(tree.getroot())
        
        return tree
        