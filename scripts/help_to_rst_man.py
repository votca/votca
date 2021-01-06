#! /usr/bin/env python3

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
import subprocess
from typing import NamedTuple


class Section(NamedTuple):
    """Data representing the help message."""
    header: str
    content: str


def parse_user_arguments() -> argparse.Namespace:
    """Read the user arguments."""
    parser = argparse.ArgumentParser("help2rstman")
    parser.add_argument("-n", "--name", required=True,
                        help="Name of the tool to extract the help message")
    parser.add_argument(
        "-f", "--format", choices=["rst", "man"], default="rst", help="Output Format"),
    parser.add_argument("-o", "--out", help="Output file")

    return parser.parse_args()


def parse_help_message(msg: str):
    """Parse the help message into an intermediate representation."""
    print(msg.splitlines())


def convert_help_message(args: argparse.Namespace):
    """Convert help message."""
    msg = subprocess.check_output(f"{args.name} --help", shell=True)
    parse_help_message(msg.decode())


def main():
    args = parse_user_arguments()
    convert_help_message(args)


if __name__ == "__main__":
    main()
