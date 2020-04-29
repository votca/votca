#!/usr/bin/env python
"""Module to extract metadata from the xml files."""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple

msg = "extract_metadata.py -i file.xml"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True,
                    help="Input file in YAML format")
parser.add_argument('-o', help="Optional output file", default=None)
parser.add_argument("-m", "--mode", help="Operation mode: xtp or csg",
                    choices=["xtp", "csg"], default="xtp")

MAXIMUM_LINE_LENGTH = 60

XTP_TABLE_HEADER = """
{}
The following table contains the input options for the calculator,

.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input""".format

CSG_TABLE_HEADER = """
{}
The following table contains the input options for CSG,

.. list-table::
   :header-rows: 1
   :align: center

   * - Property Name
     - Description
     - Default Value""".format


XTP_TABLE_LINE = """
   * - {}
     - {}
     - {}
     - {}""".format

CGS_TABLE_LINE = """
   * - {}
     - {}
     - {}""".format


def wrap_line(line: str) -> str:
    """Split a line into lines smaller than ``max_len``."""
    acc = [[]]
    count = 7  # mulitiple lines start at column number 7
    for word in line.split():
        if count > MAXIMUM_LINE_LENGTH:
            count = 7
            acc.append([])
        # Append words to the last list
        acc[-1].append(word)
        count += 1 + len(word)

    it = iter(" ".join(word for word in line) for line in acc)
    head = next(it)
    spaces = f"\n{' ':^6s} | "
    return "| " + head + spaces + f"{spaces}".join(line for line in it)


def xtp_extract_metadata(file_name: str) -> Tuple[str, List[ET.Element]]:
    """Get the description and elements from the xml file."""
    tree = ET.parse(file_name)
    root = tree.getroot()
    data = list(root)[0]
    header = data.attrib.get("help", "")
    return header, iter(data)


def csg_extract_metadata(file_name: str) -> Tuple[str, List[ET.Element]]:
    """Get the description and elements from the xml file."""
    tree = ET.parse(file_name)
    root = tree.getroot()
    children = list(root)
    header = children[0].text
    return header, children[1:]


def xtp_get_recursive_attributes(elem: ET.Element, root_name: str = "") -> str:
    """Get recursively the attributes of an ``ET.Element``."""
    s = ""
    if list(elem):
        return ''.join(xtp_get_recursive_attributes(el, f"{elem.tag}.") for el in list(elem))

    name = root_name + elem.tag
    description = elem.attrib.get("help", "")
    if len(description) > MAXIMUM_LINE_LENGTH:
        description = wrap_line(description)
    default = elem.attrib.get("default", "")
    choices = elem.attrib.get("choices", "")
    s += XTP_TABLE_LINE(name, description, default, choices)

    return s


def csg_get_recursive_attributes(elem: ET.Element, root_name: str = "") -> str:
    """Get recursively the attributes of an ``ET.Element``."""
    s = ""
    if list(elem):
        children = iter(elem)
        desc_elem = elem.find("DESC")
        if desc_elem is not None:
            description = desc_elem.text
        else:
            description = ""
        # Create multiline
        if len(description) > MAXIMUM_LINE_LENGTH:
            description = wrap_line(description)

        name = root_name + elem.tag
        default = "" if elem.text is None else ' '.join(elem.text.split())
        s += CGS_TABLE_LINE(name, description, default)
        s += ''.join(
            csg_get_recursive_attributes(el,
                                         f"{name}.") for el in children if el.tag != "DESC")
        return s

    return s


def xtp_create_rst_table(file_name: str) -> str:
    """Create an RST table using the metadata in the XML file."""
    header, elements = xtp_extract_metadata(file_name)
    s = XTP_TABLE_HEADER(header)
    for elem in elements:
        s += xtp_get_recursive_attributes(elem)
    return s


def csg_create_rst_table(file_name: str) -> str:
    """Create an RST table using the metadata in the XML file."""
    header, elements = csg_extract_metadata(file_name)
    s = CSG_TABLE_HEADER(header)
    for elem in elements:
        s += csg_get_recursive_attributes(elem)
    return s


def create_parent_folder(path: Path) -> None:
    """Create parent folder for ``path`` if it doesn't exists."""
    if not path.parent.exists():
        path.parent.mkdir(parents=True)


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    if args.mode == "xtp":
        table = xtp_create_rst_table(args.i)
    else:
        table = csg_create_rst_table(args.i)

    # Print
    if args.o is not None:
        path = Path(args.o)
        create_parent_folder(path)
        with open(path, 'w')as f:
            f.write(table)
    else:
        print(table)


if __name__ == "__main__":
    main()
