"""Module to edit the options in the XML file."""

import xml.etree.ElementTree as ET
from typing import Any, Dict


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

