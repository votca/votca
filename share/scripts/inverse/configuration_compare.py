#!/usr/bin/env python3
#
# Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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


from optparse import OptionParser

usage = "usage: %prog [options] conf1 conf2"
parser = OptionParser(usage=usage)
parser.add_option("--eps", dest="eps", metavar="EPS",
                  help="tolerance for mismatch", default=1e-2)
(options, args) = parser.parse_args()


def die(str):
    print(str)
    quit(254)


if len(args) != 2:
    die("two configurations required as parameters")

# open both files
try:
    conf1 = open(args[0])
    conf2 = open(args[1])
except:
    die("error while opening files")

# skip the first line
conf1.readline()
conf2.readline()

# compare the second line (nat)
if str(conf1.readline()).strip() != str(conf2.readline()).strip():
    die("nat does not match")

conf1 = conf1.readlines()
conf2 = conf2.readlines()
# loop through files
for i, (c1, c2) in enumerate(zip(conf1, conf2)):
    vals1 = c1.split()
    vals2 = c2.split()
    # 6 cols expected
    if len(vals1) != 6:
        break
    # compare 1st, 2nd, 3rd column without tolerance
    for j in range(3):
        if vals1[j] != vals2[j]:
            die("mismatch in line "+str(i+3)+" col " +
                str(j+1)+", "+str(vals1[j])+"!="+str(vals2[j]))

    # compare 4th-6th col with tolerance
    for j in range(3, 6):
        if abs(float(vals1[j]) - float(vals2[j])) > options.eps:
            die("mismatch in line "+str(i+3)+" col " +
                str(j+1)+", "+str(vals1[j])+"!="+str(vals2[j]))

# compare last line without tolerance
for j in range(3):
    if vals1[j] != vals2[j]:
        die("mismatch in line "+str(i+3)+" col " +
            str(j+1)+", "+str(vals1[j])+"!="+str(vals2[j]))
