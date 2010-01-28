#!/usr/bin/env python

import sys

import pygsl.sf
import getopt
from pygsl import spline
from pygsl import _numobj as numx
import math

xvalues = []
yvalues = []



outfile = ""

doresample =False
dosmoothtorho_0 = False
dowritedpotf = False

def weight(x):
    c = math.cos(math.pi/(2*(xstop-xstart)*0.1)*x)
    return c*c


options = ["xstart=", "xstop=", "step=", "infile=", "outfile="]

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


        
        
myspline = spline.cspline(len(xvalues))
myspline.init(xvalues,yvalues)
f = open(outfile,"w")


i = 0
tempx = []
tempy = []
for x in xvalues:
    tempx.append (x)
    if x-xstart < 0.1*(xstop-xstart): 

        tempy.append ((1-weight(math.fabs(x-xstart)))*yvalues[i])
    elif x-xstart > 0.9*(xstop-xstart): 

        tempy.append ((1-weight(math.fabs(xstop-x)))*yvalues[i])
    else:
        tempy.append(yvalues[i])
    i=i+1

i = 0
for x in tempx:
    f.write('%15.10e %15.10e i\n' % (x-xstart, tempy[i]))
    i=i+1
