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
import numpy
from read_orbitals import *

#path_scforb_1=""
#path_scforb_2=""
path_molA=""
path_molB=""



arguments = sys.argv

if len(arguments) != 3:
 print 'merge_orbitals.py called with insuffcient number of arguments'
 sys.exit()

#specification of monomer 1
path_molA=arguments[1]
scfconv,mono1Dim=getOrbDim(path_molA+"/mos")
scfconv1=scfconv

#specification of monomer 2
path_molB=arguments[2]
scfconv,mono2Dim=getOrbDim(path_molB+"/mos")
scfconv2=scfconv

#read orbitals for monomer 1
if path_molA!="":
 path=path_molA+'/mos'
 vec_mos_1=readOrb(path,1,mono2Dim)
 
 
#read orbitals for monomer 2
if path_molB!="":
 path=path_molB+'/mos'
 vec_mos_2=readOrb(path,2,mono1Dim)

nsaosA=len(vec_mos_1)
nsaosB=len(vec_mos_2)
nsaosD=int(nsaosA)+int(nsaosB)

#now merge mos
vec_mos_1.extend(vec_mos_2)
merged_mos=vec_mos_1

#do sorting and printout
merged_mos.sort()
nsaos=len(merged_mos)
writeMo(scfconv,nsaos,merged_mos,"mos")


sys.exit()






