#!/usr/bin/env python
"""Module to extract metadata from the xml files."""

import argparse
import xml.etree.ElementTree as ET
from typing import Generator, List

msg = "extract_metadata.py -i file.xml"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True,
                    help="Input file in YAML format")

TABLE_HEADER = """
+------------------------------+------------------------------------+-------------------+--------------------+
|  Property Name               |  Description                       | Default Value     |   Valid Input      |
+==============================+====================================+===================+====================+
"""

TABLE_SEPARATOR = """+------------------------------+------------------------------------+-------------------+--------------------+\n"""

LEN_DESCRIPTION = 36
TABLE_LINE = """|{:^30}|{:^36}|{:^19}|                    |\n""".format


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
    if len(description) > LEN_DESCRIPTION:
        s += wrap_line(name, description, default)
    else:
        s += TABLE_LINE(name, description, default)
    s += TABLE_SEPARATOR

    return s


def create_rst_table(file_name: str) -> None:
    """Create an RST table using the metadata in the XML file."""
    s = TABLE_HEADER
    for elem in extract_metadata(file_name):
        s += get_recursive_attributes(elem)
    print(s)


def wrap_line(name: str, description: str, default: str) -> str:
    """Split the line in several ones."""
    xs = split_lines(description)
    s = TABLE_LINE(name, xs[0], default)
    for x in xs[1:]:
        # | and 30 white spaces follow by |`x`| follow by | and 19 white spaces and |20 spaces|
        w = f"|{30 * ' '}|{x:<36}|{19 * ' '}|{20 * ' '}|\n"
        s += w

    return s


def split_lines(data: str) -> List[str]:
    """Split the ``data`` in lines smaller that ``LEN_DESCRIPTION``.

    It Also appends the | char at the beginning to force the rst
    interpreter to split the line.
    """
    acc = [[]]
    # there are 3 characters in " | "
    count = 3
    for x in data.split():
        # create new line
        if count + len(x) + 1 > LEN_DESCRIPTION:
            count = 3
            acc.append([])
        acc[-1].append(x)
        count += len(x) + 1

    return [" | " + ' '.join(x) for x in acc]


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    create_rst_table(args.i)


if __name__ == "__main__":
    main()
