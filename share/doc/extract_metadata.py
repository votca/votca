#!/usr/bin/env python
"""Module to extract metadata from the xml files."""

import argparse
import xml.etree.ElementTree as ET
from typing import Generator

msg = "extract_metadata.py -i file.xml"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True,
                    help="Input file in YAML format")
parser.add_argument('-o', help="Optional output file", default=None)

TABLE_HEADER = """
.. list-table:: Description
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input"""

TABLE_LINE = """
   * - {}
     - {}
     - {}
     - {}""".format


def extract_metadata(file_name: str) -> Generator[ET.Element, None, None]:
    """Get the metadata from the xml."""
    tree = ET.parse(file_name)
    root = tree.getroot()
    data = root.getchildren()[0]
    for elem in data.getchildren():
        yield elem


def get_recursive_attributes(elem: ET.Element, root_name: str = "") -> str:
    """Get recursively the attributes of an ``ET.Element``."""
    s = ""
    if elem.getchildren():
        return ''.join(get_recursive_attributes(el, f"{elem.tag}.") for el in elem.getchildren())

    name = root_name + elem.tag
    description = elem.attrib.get("help", "")
    default = elem.attrib.get("default", "")
    choices = elem.attrib.get("choices", "")
    s += TABLE_LINE(name, description, default, choices)

    return s


def create_rst_table(file_name: str) -> str:
    """Create an RST table using the metadata in the XML file."""
    s = TABLE_HEADER
    for elem in extract_metadata(file_name):
        s += get_recursive_attributes(elem)
    return s


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    table = create_rst_table(args.i)
    if args.o is not None:
        with open(args.o, 'w')as f:
            f.write(table)
    else:
        print(table)


if __name__ == "__main__":
    main()
