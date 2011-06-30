#!/usr/bin/env python

import sys
import numpy
from read_orbitals import *

path_uhforb_1=""
path_scforb_1=""
path_uhforb_2=""
path_scforb_2=""

arguments = sys.argv

if len(arguments) != 3:
 print 'usage: merge_own.py <directory of monomer 1> <directory of monomer 2>'
 sys.exit()

#specification of monomer 1
path_scforb_1=arguments[1]
scfconv,mono1Dim=getOrbDim(path_scforb_1+"/mos")
scfconv1=scfconv

#specification of monomer 2
path_scforb_2=arguments[2]
scfconv,mono2Dim=getOrbDim(path_scforb_2+"/mos")
scfconv2=scfconv

#read orbitals for monomer 1
if path_scforb_1!="":
 path=path_scforb_1+'/mos'
 vec_mos_1=readOrb(path,1,mono2Dim)
 
 
#read orbitals for monomer 2
if path_scforb_2!="":
 path=path_scforb_2+'/mos'
 vec_mos_2=readOrb(path,2,mono1Dim)

nsaosA=len(vec_mos_1)
nsaosB=len(vec_mos_2)
nsaosD=int(nsaosA)+int(nsaosB)

#now merge mos
vec_mos_1.extend(vec_mos_2)
merged_mos=vec_mos_1

#do sorting and printout
merged_mos.sort()
path=path_scforb_1+'/mos'
nsaos=len(merged_mos)
writeMo(scfconv,nsaos,merged_mos,"mos")


sys.exit()






