#!/usr/bin/env python2
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
import os

#import pygsl.sf
import getopt
#from pygsl import spline
#from pygsl import _numobj as numx
import math

xvalues = []
yvalues = []

outfile = ""

options = ["xstart=", "xstop=", "infile=", "outfile=","help"]

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
This script get the a subset of a table

Usage: %(name)s 
Allowed options:
    --xstart     X.X  x value where the subset starts
    --xstop      X.X  x value where the subset stops
    --infile    FILE  input file
    --outfile   FILE  output file
""" % {'name': os.path.basename(sys.argv[0]),'ver': '%version%'}
      sys.exit(0)
    elif o == "-v":
        verbose = True
    elif o == "--xstart":
        xstart = float(a)
    elif o == "--xstop":
        xstop = float(a)
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
		if float(values[0]) >= xstart and float(values[0]) <= xstop:
			xvalues.append(float(values[0]))
			yvalues.append(float(values[1]))

f = open(outfile,"w")


i = 0
tempx = []
tempy = []
for x in xvalues:
    tempx.append (x)
    tempy.append(yvalues[i])
    i=i+1

i = 0
for x in tempx:
    f.write('%15.10e %15.10e i\n' % (x-xstart, tempy[i]))
    i=i+1
