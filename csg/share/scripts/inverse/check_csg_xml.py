#!/usr/bin/env python3
"""Script to check csg XML file for invalid tags."""
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

import argparse
import sys
import xml.etree.ElementTree as ET

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")


def main():
    # get command line arguments
    args = get_args()
    # check XML for invalid tags
    check_votca_settings_xml(
        ET.fromstring(args.csg_xml_file.read()),
        ET.fromstring(args.csg_xml_defaults_file.read()),
    )


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Check the csg XML file for invalid tags.

    It works by comparing to the tags in the default XML file.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "csg_xml_file",
        type=argparse.FileType("r"),
        help="csg XML file",
    )
    parser.add_argument(
        "csg_xml_defaults_file",
        type=argparse.FileType("r"),
        help="csg XML defaults file",
    )
    # parse
    args = parser.parse_args()
    return args


def check_votca_settings_xml(root, root_defaults):
    def iter_xml(node, path="."):
        yield node, path
        for child in node:
            child_path = f"{path}/{child.tag}"
            for child, child_path in iter_xml(child, path=child_path):
                yield child, child_path

    found_bad_tags = []
    for child, child_path in iter_xml(root):
        # options for cg.bonded are partially listed in cg.non-bonded
        child_path_non_bonded = child_path.replace("/bonded/", "/non-bonded/")
        if (
            root_defaults.find(child_path) is None
            and root_defaults.find(child_path_non_bonded) is None
            # per group settings can have any name
            and not child_path.startswith("./cg/inverse/imc/")
        ):
            found_bad_tags.append(child_path)
    if len(found_bad_tags) > 0:
        found_bad_tags_lines = "\n".join(found_bad_tags)
        print(
            f"The settings XML file contains the tags\n\n{found_bad_tags_lines}\n\nbut "
            "those paths does not exist in the XML defaults file and are therefore "
            "not supported."
        )


if __name__ == "__main__":
    main()
