#!/usr/bin/env python
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

import sys
import getopt
import math

firstx = []
firsty = []
secondx = []
secondy = []

outfile = ""

options = ["adressc=", "infile=", "outfile=","help"]

try:
    opts, args = getopt.getopt(sys.argv[1:], "", options)
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    print options
    sys.exit(2)
for o, a in opts:
    if o == "--help":
      print """%(name)s, version %(ver)s 
This script symmetrizes the density around --adressc for thermodynamic force iteration

Usage: %(name)s 
Allowed options:
    --adressc    X.X  center of the adress zone (x-value)
    --infile    FILE  input file
    --outfile   FILE  output file
""" % {'name': sys.argv[0],'ver': '%version%'}
      sys.exit(0)
    elif o == "-v":
        verbose = True
    elif o == "--adressc":
        adressc = float(a)
    elif o in ("--infile"):
        infile = a
    elif o in ("--outfile"):
        outfile = a
    else:
        print options
        assert False, "unhandled option"

for line in open(infile,"r").readlines():
	if line[0] != "@" and line[0] != "#":
		values = line.split()
		if float(values[0]) <= adressc:
			firstx.append(float(values[0]))
			firsty.append(float(values[1]))
                        if float(values[0]) == adressc:
                            secondx.append(float(values[0]))
                            secondy.append(float(values[1]))
                else:
                       
                        if len(firstx)-1-len(secondx) >= 0 and (len(firstx)-1-len(secondx)) < len (firsty):
                            secondx.append(float(values[0])) 
			    secondy.append( 0.5*(firsty[len(firstx)-len(secondx)]+(float(values[1]))) )
                        else:
                            print "Warning: symmetrize_density.pl : adressc not in center of data", line
                            print "index", len(firstx)-len(secondx)

f = open(outfile,"w")
i=0
for x in secondx:
    f.write('%15.10e %15.10e i\n' % (x, secondy[i]))
    i=i+1
