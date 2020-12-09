#!/usr/bin/env python
#
# Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Module to extract metadata from the xml files."""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

msg = "extract_metadata.py -i file.xml"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True,
                    help="Input file in YAML format")
parser.add_argument('-o', help="Optional output file", default=None)
parser.add_argument("-m", "--mode", help="Operation mode: xtp, csg, qm",
                    choices=["xtp", "csg", "qm"], default="xtp")

MAXIMUM_LINE_LENGTH = 60

DEFAULTS_TABLE_HEADER = """
.. list-table::
   :header-rows: 1
   :widths: 30 20 15 15
   :align: center

   * - Property Name
     - Description
     - Default Value
     - Valid Input"""


def xtp_table_header(x: str) -> str:
    """Create a table for XTP calculators."""
    return f"""
{x}
The following table contains the defaults input options for the calculator,
{DEFAULTS_TABLE_HEADER}"""


def csg_table_header(x: str):
    """Create a CSG table with its default values."""
    return f"""
{x}
The following table contains the input options for CSG,

.. list-table::
   :header-rows: 1
   :align: center

   * - Property Name
     - Description
     - Default Value"""


XTP_TABLE_LINE = """
   * - {}
     - {}
     - {}
     - {}""".format

CGS_TABLE_LINE = """
   * - {}
     - {}
     - {}""".format


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    file_name = Path(args.i)
    if args.mode == "xtp":
        table = xtp_create_rst_table(file_name)
    elif args.mode == "csg":
        table = csg_create_rst_table(file_name)
    else:
        table = qmpackage_create_rst_table(file_name)

    # Print
    if args.o is not None:
        path = Path(args.o)
        create_parent_folder(path)
        with open(path, 'w')as f:
            f.write(table)
    else:
        print(table)


def get_root_children(file_name: Path) -> List[ET.Element]:
    """Get all the node children from root."""
    tree = ET.parse(file_name)
    root = tree.getroot()
    return list(root)


def xtp_extract_metadata(file_name: Path) -> Tuple[str, ET.Element]:
    """Get the description and elements from the xml file."""
    data = get_root_children(file_name)[0]
    header = data.attrib.get("help", "")
    return header, data


def csg_extract_metadata(file_name: Path) -> Tuple[str, List[ET.Element]]:
    """Get the description and elements from the xml file."""
    children = get_root_children(file_name)
    header = children[0].get("text", "")
    return header, children[1:]


def xtp_get_recursive_attributes(elem: ET.Element, root_name: str = "") -> str:
    """Get recursively the attributes of an ``ET.Element``."""
    s = ""
    name = root_name + elem.tag
    if list(elem):
        return ''.join(xtp_get_recursive_attributes(el, f"{name}.") for el in list(elem))

    description = split_line(elem.attrib.get("help", ""))
    default = split_line(elem.attrib.get("default", ""))
    choices = multiline(elem.attrib.get("choices", ""))
    s += XTP_TABLE_LINE(name, description, default, choices)

    return s


def csg_get_recursive_attributes(elem: ET.Element, root_name: str = "") -> str:
    """Get recursively the attributes of an ``ET.Element``."""
    s = ""
    if list(elem):
        children = iter(elem)
        desc_elem = elem.find("DESC")
        if desc_elem is not None:
            description = desc_elem.get("text", "")
        else:
            description = ""
        description = split_line(description)

        name = root_name + elem.tag
        default = "" if elem.text is None else ' '.join(elem.text.split())
        s += CGS_TABLE_LINE(name, description, default)
        s += ''.join(
            csg_get_recursive_attributes(el,
                                         f"{name}.") for el in children if el.tag != "DESC")
        return s

    return s


def xtp_create_rst_table(file_name: Path) -> str:
    """Create an RST table using the metadata in the XML file."""
    header, elements = xtp_extract_metadata(file_name)
    header = generate_title(file_name.stem) + header
    s = xtp_table_header(header) if elements else f"{header}\n"
    for elem in elements:
        s += xtp_get_recursive_attributes(elem)

    s += generate_note(file_name.stem)
    return s


def csg_create_rst_table(file_name: Path) -> str:
    """Create an RST table using the metadata in the XML file."""
    header, elements = csg_extract_metadata(file_name)
    s = csg_table_header(header)
    for elem in elements:
        s += csg_get_recursive_attributes(elem)
    return s


def qmpackage_create_rst_table(file_name) -> str:
    """Create an RST table using the metadata in the XML file."""
    children = get_root_children(file_name)
    s = DEFAULTS_TABLE_HEADER
    for elem in children:
        s += xtp_get_recursive_attributes(elem)
    return s


def create_parent_folder(path: Path) -> None:
    """Create parent folder for ``path`` if it doesn't exists."""
    if not path.parent.exists():
        path.parent.mkdir(parents=True)


def split_line(line: str, sep: Optional[str] = None) -> str:
    """Split line if larger than ``MAXIMUM_LINE_LENGTH``."""
    if len(line) > MAXIMUM_LINE_LENGTH:
        return wrap_line(line, sep)
    else:
        return line


def wrap_line(line: str, sep: Optional[str]) -> str:
    """Split a line into lines smaller than ``max_len``."""
    acc = [[]]  # type: List[List[str]]
    count = 7  # mulitiple lines start at column number 7
    for word in line.split(sep=sep):
        # If cumulative sum or the length of word is greater than MAXIMUM_LINE_LENGTH
        if any(x > MAXIMUM_LINE_LENGTH for x in (count, len(word))):
            count = 7
            acc.append([])
        # Append words to the last list
        acc[-1].append(word)
        count += 1 + len(word)

    it = iter(" ".join(word for word in line) for line in acc)
    return column_multiline(it)


def generate_title(stem: str) -> str:
    """Generate title in rst format using ``file_name``."""
    return f"{stem}\n{'*' * len(stem)}\n"


def generate_note(stem: str) -> str:
    """Generate note specifying path to the xml file."""
    note = f"""
.. note::
   An *xml* file containing the defaults for the `{stem}` calculator can be found at `$VOTCASHARE/xtp/xml/{stem}.xml`
"""
    return note


def multiline(line: str) -> str:
    """Split the comma separated words into lines."""
    xs = line.split(',')
    if len(xs) >= 2:
        return column_multiline(iter(xs))
    else:
        return line


def column_multiline(xs: Iterable[str]) -> str:
    """Create a column with multiple lines."""
    spaces = f"\n{' ':^6s} | "
    return "| " + next(xs) + spaces + f"{spaces}".join(xs)


if __name__ == "__main__":
    main()
