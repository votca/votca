"""Module defining the Option class."""
from typing import Any, Dict
import os
import xml.etree.ElementTree as ET
from lxml import etree
from lxml.etree import Element
from pyxtp.xml_editor import NestedNamespace, xml2namespace, namespace2dict

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

        votcashare = os.environ.get('VOTCASHARE')
        if votcashare is None:
            raise ValueError('Please set VOTCASHARE environement variable')

        self.dftgwbse_options_file = f'{votcashare}/xtp/xml/dftgwbse.xml'
        self.gwbse_options_file = f'{votcashare}/xtp/xml/gwbse.xml'
        self.dftpackage_options_file = f'{votcashare}/xtp/xml/dftpackage.xml'

        self.data = xml2namespace(self.dftgwbse_options_file).options.dftgwbse
        self.data.dftpackage = xml2namespace(self.dftpackage_options_file).dftpackage
        self.data.gwbse = xml2namespace(self.gwbse_options_file).gwbse
      
        
    def write_xml(self):
        """Writes the XML files containing the user defined options
        """

        # create a dict from the user defined options
        options = namespace2dict(self.data)
        
        # creates sub dict specific to the different xml files
        dftgwbse_options = {k:v for k,v in options.items() if not any(k.startswith(s) for s in ['dftpackage','gwbse'])}
        dft_options = {k.lstrip('dftpackage/'):v for k,v in options.items() if k.startswith('dftpackage')}
        gwbse_options = {k.lstrip('gwbse/'):v for k,v in options.items() if k.startswith('gwbse')}
        
        # update the xml files with the user defined options
        xml_dftgwbse = self._update(self.dftgwbse_options_file, dftgwbse_options)
        xml_dft = self._update(self.dftpackage_options_file, dft_options)
        xml_gwbse = self._update(self.gwbse_options_file, gwbse_options)
        
        # write the xml files
        xml_dftgwbse.write('dftgwbse.xml')
        xml_dft.write('dftpackage.xml')
        xml_gwbse.write('gwbse.xml')
        
        # clean/write the data
        self._clean('dftgwbse.xml')
        self._clean('dftpackage.xml')
        self._clean('gwbse.xml')

    @staticmethod
    def _update(xml_filename: str, dict_options: dict) -> ET.ElementTree:
        """Read and update default option files with the user define options

        Args:
            xml_filename (str): default XML file
            dict_options (dict): dictionary of user defined options

        Returns:
            ET.ElementTree: xml structure containing the update options.
        """
        # read the xml file
        xml = ET.parse(xml_filename)
        
        # replace text value of the dict elements
        for key, value in dict_options.items():
            try:
                xml.find('.//'+key).text = value
            except:
                print(key, ' not in ', xml_filename)
        return xml 
    
    @staticmethod
    def _clean(xml_filename: str, remove_attributes=[]):
        """Recursively remove all empty nodes and elements containing only empty nodes

        Args:
            xml_filename (str): name of the xml file
            remove_attributes (list): list of xml attributes to remove
        """
        
        def is_empty(elem: Element) -> bool:
            """returns true if e does not have text

            Args:
                elem (Element): xml element

            Returns:
                bool true if e is empty
            """
            return (elem.text is None) or (elem.text.strip() == '')
        
        def recursively_remove_empty(elem: Element):
            """_summary_

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
            
        # parse the file
        tree = etree.parse(xml_filename)
        # remove attributes if any are passed
        etree.strip_attributes(tree, remove_attributes)
        # remove the empties
        recursively_remove_empty(tree.getroot())
        # write to file
        tree.write(xml_filename)