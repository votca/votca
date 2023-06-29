"""Module to edit the options in the XML file."""

import xml.etree.ElementTree as ET
import os
from lxml import etree
from lxml.etree import Element, ElementTree
from copy import deepcopy

def update_xml_text(xml: ElementTree, dict_options: dict) -> ElementTree:
    """update default option files with the user define options

    Args:
        xml_tree (ElementTree): tree of the xml file
        dict_options (dict): dictionary of user defined options

    Returns:
        ET.ElementTree: xml structure containing the update options.
    """

    # copy the xml tree
    xml_cpy = deepcopy(xml)

    # replace text value of the dict elements
    for key, value in dict_options.items():
        try:
            xml_cpy.find('.//'+key).text = str(value)
        except Exception as e:
            print(e)
            print(key, ' not found ')
    return xml_cpy


def insert_linked_xmlfiles(tree: ElementTree, base_path: str) -> None:
    """insert the tree of linked xml files in the original tree

    Args:
        tree (ElementTree): original xml tree containing linked files
        base_path (str): base path for the linked files
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
            link_name = os.path.join('subpackages', el.attrib['link'])
            link_name = os.path.join(base_path, link_name)
            if not os.path.isfile(link_name):
                raise ValueError('File %s not found' %link_name)
            xml_link = etree.parse(link_name)
            el.getparent().replace(el, xml_link.getroot())

    # recursively insert the links
    _recursive_insert_links(tree.getroot())

    return tree


def remove_empty_text_elements(tree: ET.ElementTree, remove_attributes=[]) -> ElementTree:
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

    def recursively_remove_empty_elements(elem: Element):
        """reccursive function that remove empty elements

        Args:
            e (Element): xml element

        Returns:
            bool: True if we need to remove the element False otherwise
        """

        # recursive call that gets a list of :
        # [] if e doens't have children
        # [True, False, False, True ] where each bool indicate if e children are empty or not
        child_bool = list(recursively_remove_empty_elements(c) for c in elem.iterchildren())
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
    recursively_remove_empty_elements(tree.getroot())

    return tree