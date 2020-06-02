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

"""Module to edit inplace the content of xml files."""

__all__ = ["edit_calculator", "add_section"]

import io
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple


def get_tree_and_root(file_name: Path) -> Tuple[ET.Element, ET.Element]:
    """Get the tree and root of a calculator input."""
    if not file_name.exists():
        raise RuntimeError(f"There is no file: {file_name}")
    # Read XML file
    tree = ET.parse(file_name)
    root = tree.getroot()

    return tree, root


def edit_calculator(calculator: str, option: str, value: str, folder: str = "OPTIONFILES") -> None:
    """Edit ``value`` from ``option`` of calculator."""
    file_name = Path(folder) / Path(f"{calculator}.xml")
    tree, root = get_tree_and_root(file_name)
    edit_option(list(root)[0], option, value)
    # Finally overwrite the file
    tree.write(file_name)
    msg = f"The option '{option}' on file '{file_name.as_posix()}' has been set to '{value}'"
    print(msg)


def edit_option(elem: ET.Element, option: str, value: str) -> None:
    """Change value of ``option`` in ``elem``."""
    children = find_nodes(elem, option)
    if not children:
        raise RuntimeError(f"There is not {option} in {elem}")
    else:
        for child in children:
            child.text = value


def find_nodes(elem: ET.Element, option: str) -> List[ET.Element]:
    """Find ``option`` in ``elem`` and in all its children."""
    sections = option.split('.')
    if len(sections) == 1:
        acc = []
        find_nodes_recursively(elem, option, acc)
        return acc
    else:
        return find_section(elem, sections)


def find_section(elem: ET.Element, sections: List[str]) -> List[ET.Element]:
    """Search recursively for '.' separated section in elem."""
    for name in sections:
        elem = elem.find(name)
        if elem is None:
            break

    return [elem]


def find_nodes_recursively(elem: ET.Element, option: str, acc: List[ET.Element]) -> None:
    """Find ``option`` in ``elem`` and in all its children."""
    child = elem.find(option)
    if child is not None:
        acc.append(child)
    # If elem has children
    if list(elem):
        for c in list(elem):
            find_nodes_recursively(c, option, acc)


def add_section(calculator: str, xml_string: str, folder: str = "OPTIONFILES") -> None:
    """Add a new element to the ``calculator`` input."""
    file_name = Path(folder) / Path(f"{calculator}.xml")
    tree, root = get_tree_and_root(file_name)
    new_elem = ET.parse(io.StringIO(xml_string))
    calculator = next(iter(root))
    calculator.append(new_elem.getroot())
    tree.write(file_name)
