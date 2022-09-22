"""Module defining the Option class."""
from typing import Any, Dict
import os
import xml.etree.ElementTree as ET
from lxml import etree
from lxml.etree import Element, ElementTree

from pyxtp.nestednamespace import NestedNamespace
from pyxtp.xml_editor import remove_empty_text_elements, insert_linked_xmlfiles, update_xml_text


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

        # main filename
        dftgwbse_options_file = f'{votcashare}/xtp/xml/dftgwbse.xml'      
        
        # parse xmlfile and replace links to the embedded xml files
        self._xml_tree = insert_linked_xmlfiles(etree.parse(dftgwbse_options_file), base_path=f'{votcashare}/xtp/xml/')
        
        # create a namespace from the xml tree
        self.__dict__.update(NestedNamespace._fromxml(etree.tostring(self._xml_tree)).options.dftgwbse.__dict__)
    

    
    def _write_xml(self, filename: str) -> None:
        """update and clean the xml and then write the input file

        Args:
            filename (str): filename in the 
        """
        
        # create a dict of the user providedoptions
        dftgwbse_options = self._todict()
        
        # update the options from the user provided input
        self._xml_tree = update_xml_text(self._xml_tree, dftgwbse_options)

        # remove the empty elements
        self._xml_tree = remove_empty_text_elements(self._xml_tree)
        
        # write the file
        self._xml_tree.write(filename)
    

        