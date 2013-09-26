#!/usr/bin/env python2
## his2xyz.py
## function: convert the HISTORY file to the multi-frame XYZ file
## Jen-Chang Chen
## version 1.0
## Jan 25, 2005
## Usage: 
## I. chmod 755 his2xyz.py 
## II. ./his2xyz.py HISTORY HIS.xyz

import sys,string

atomList=['A','B','CG']
inputFile=open(sys.argv[1],'r')
outFile=open(sys.argv[2],'w')

title=inputFile.readline()
line=inputFile.readline()

while(line!=""): 
	if string.split(line)[0]=='timestep':
		timestep='step= '+string.split(line)[1]+'\n'
		totalAtom=string.split(line)[2]+'\n'	
		outFile.write(totalAtom)
		outFile.write(timestep)
		
	if string.split(line)[0][:2] in atomList:
			atomName=string.split(line)[0][:2]+'\t'
			outFile.write(atomName)
			xyzline=inputFile.readline()
			outFile.write(xyzline)

	line=inputFile.readline()

inputFile.close()
outFile.close()
