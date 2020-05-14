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

__all__ = ["edit_calculator"]

import io
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Tuple


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
    child = find_node_recursively(elem, option.lower())
    child.text = value


def find_node_recursively(elem: ET.Element, option: str) -> ET.Element:
    """Check if ``elem`` of its children has a node called ``option``."""
    child = elem.find(option.lower())
    if child is not None:
        return child
    elif list(elem):
        for c in list(elem):
            x = find_node_recursively(c, option)
            if x is not None:
                return x
    else:
        return None


def add_section(calculator: str, xml_string: str, folder: str = "OPTIONFILES") -> None:
    """Add a new element to the ``calculator`` input."""
    file_name = Path(folder) / Path(f"{calculator}.xml")
    tree, root = get_tree_and_root(file_name)
    new_elem = ET.parse(io.StringIO(xml_string))
    root.append(new_elem.getroot())
    tree.write(file_name)
