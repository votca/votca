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

if not sys.version_info >= (3, 6):
    raise Exception("This script needs Python 3.6+.")

REPORT_HEADER = (
    "The settings XML file contains the tags which do not exist in the XML defaults "
    "file:"
)


def main():
    # get command line arguments
    args = get_args()
    # check XML for invalid tags
    report = check_votca_settings_xml(
        args.csg_xml_file.read(),
        args.csg_xml_defaults_file.read(),
    )
    print(report)


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


def check_votca_settings_xml(root_xml, root_defaults_xml):
    """Check xml tags in root against root_defaults"""

    # converte to ElementTree
    root = ET.fromstring(root_xml)
    root_defaults = ET.fromstring(root_defaults_xml)

    def iter_xml(node, path="."):
        """Recurse an ElementTree"""
        yield node, path
        for child in node:
            child_path = f"{path}/{child.tag}"
            for child, child_path in iter_xml(child, path=child_path):
                yield child, child_path

    sim_progs = [
        "gromacs",
        "lammps",
        "hoomd",
        "espresso",
        "dlpoly",
        "espressopp",
        "hoomd-blue",
    ]
    bad_paths = []
    for _, child_path in iter_xml(root):
        # options for cg.bonded are partially listed in cg.non-bonded
        child_path_non_bonded = child_path.replace("/bonded/", "/non-bonded/")
        # options for cg.$program are listed in cg.sim_prog, similar per interaction
        child_path_sim_prog = child_path
        for sim_prog in sim_progs:
            child_path_sim_prog = child_path_sim_prog.replace(
                f"inverse/{sim_prog}/", "inverse/sim_prog/"
            )
        if (
            root_defaults.find(child_path) is None
            and root_defaults.find(child_path_non_bonded) is None
            and root_defaults.find(child_path_sim_prog) is None
            # per group settings can have any name
            and not child_path.startswith("./inverse/imc/")
            # multidir gromacs allows multiple index, mdp, topol_in, or conf
            and not any(
                child_path.startswith(f"./inverse/gromacs/{tag}/")
                for tag in ["index", "mdp", "topol_in", "conf"]
            )
        ):
            bad_paths.append(child_path)

    # find what the user maybe meant
    suggestions = []
    for found_bad_path in bad_paths:
        last_tag = found_bad_path.split("/")[-1]
        tag_suggestions = [
            default_path.replace("/non-bonded", "/{bonded,non-bonded}")
            for default_element, default_path in iter_xml(root_defaults)
            if default_element.tag == last_tag
        ]
        suggestions.append(tag_suggestions)

    # output bad tags and suggestions
    if len(bad_paths) > 0:
        found_bad_lines = ""  # lines in output about bad paths
        for bad_path, tag_suggestions in zip(bad_paths, suggestions):
            found_bad_lines += bad_path
            if len(tag_suggestions) > 0:
                found_bad_lines += f"  ( did you mean {' or '.join(tag_suggestions)} )"
            found_bad_lines += "\n"
        return f"{REPORT_HEADER}\n{found_bad_lines}"
    else:
        return ""


def test_check_votca_settings_xml():
    report = check_votca_settings_xml("<cg></cg>", "<cg></cg>")
    assert report == ""
    report = check_votca_settings_xml("<cg><foo></foo></cg>", "<cg></cg>")
    assert report == f"{REPORT_HEADER}\n./foo\n"
    report = check_votca_settings_xml(
        "<cg><foo></foo></cg>", "<cg><bar><foo></foo></bar></cg>"
    )
    assert report == f"{REPORT_HEADER}\n./foo  ( did you mean ./bar/foo )\n"


if __name__ == "__main__":
    main()
