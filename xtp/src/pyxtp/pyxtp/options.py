"""Module defining the Option class."""
import os
from lxml import etree
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
        
        # insert the namesapce data into the instance
        self.__dict__.update(NestedNamespace.__fromxml__(etree.tostring(self._xml_tree)).options.dftgwbse.__dict__)
    
    
    def _write_xml(self, filename: str) -> None:
        """update and clean the xml and then write the input file

        Args:
            filename (str): filename in the 
        """
        
        # create a dict of the user providedoptions
        dftgwbse_options = self.__todict__()
        
        # update the options from the user provided input
        self._xml_tree = update_xml_text(self._xml_tree, dftgwbse_options)

        # remove the empty elements
        self._xml_tree = remove_empty_text_elements(self._xml_tree)
        
        # write the file
        self._xml_tree.write(filename)
        
    def _fill_default_values(self) -> None:
        
        """Fill the default values of the nested namespace.
        
        Returns:
            dict: dictionary maps the namespace
        """
        
        def _recursive_filldefault(namespace_data: NestedNamespace) -> dict:
            """Fill the default values of the nested namespace.

            Args:
                namespace_data (NestedNamespace): data in namespace form
                path (str, optional): basepath in the structure. Defaults to ''.
            
            Returns:
                dict: dictionary maps the namespace
            """

            default_dict = dict()
            for key, value in namespace_data.__dict__.items():
                
                if isinstance(value, NestedNamespace):

                    if value._haschildren():
                        _recursive_filldefault(value)
                        
                    elif '_default' in value.__dict__.keys():
                        namespace_data.__dict__[key] = value.__dict__['_default']
                        

            
        _recursive_filldefault(self)
    

        