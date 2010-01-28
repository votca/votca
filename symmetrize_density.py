#!/usr/bin/env python

import sys

import pygsl.sf
import getopt
from pygsl import spline
from pygsl import _numobj as numx
import math

firstx = []
firsty = []
secondx = []
secondy = []

outfile = ""

options = ["adressc=", "infile=", "outfile="]

try:
    opts, args = getopt.getopt(sys.argv[1:], "", options)
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    print options
    sys.exit(2)
for o, a in opts:
    if o == "-v":
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
		if float(values[0]) < adressc:
			firstx.append(float(values[0]))
			firsty.append(float(values[1]))
                else:
                        secondx.append(float(values[0]))
                        if len(firstx)-len(secondx) >= 0 and (len(firstx)-len(secondx)) < len (firsty):
                            secondy.append( 0.5*(firsty[len(firstx)-len(secondx)]+(float(values[1]))) )
                        else:
                            print "Warning: symmetrize_density.pl : adressc not in center of data", line
                            print "index", len(firstx)-len(secondx)

f = open(outfile,"w")
i=0
for x in secondx:
    f.write('%15.10e %15.10e i\n' % (x, secondy[i]))
    i=i+1
