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


#option parsing
#TODO optparse is depreceated since 2.7
#switch to argparse (NOT supported in 2.6)

from optparse import OptionParser
import math #needed for NaN checks

usage = "usage: %prog [options] conf1 conf2"
parser = OptionParser(usage=usage)
parser.add_option("--eps", dest="eps", metavar="EPS",
                  help="tolerance for mismatch", default=1e-5)
(options, args) = parser.parse_args()                  

eps = float(options.eps)

def die(str):
  print str
  quit(1)


def difference_relative(v1, v2, _eps):
  """calculates relative difference - see table_combine.pl for original version"""
  """this is just a translation of the perl function"""
  v1 = float(v1)
  v2 = float(v2)
  if math.isnan(v1): v1 = 0
  if math.isnan(v2): v2 = 0
  if v1 == v2: return 0
  if abs(v1-v2) < _eps: return abs(v1-v2)
  if abs(v1) > abs(v2): return abs(v1-v2)/abs(v1)
  return abs(v1-v2)/abs(v2)

if len(args) != 2:
  die("two data files required as parameters")

#open both files
try: 
  dataf1 = open(args[0])
  dataf2 = open(args[1])
except:
  die("error while opening files")

data1 = dataf1.readlines()
data2 = dataf2.readlines()
#loop through files
for i in range(len(data1)):
  #TODO: skip comment lines?
  #if ^[@#]: continue

  vals1 = data1[i].split()
  vals2 = data2[i].split()
  
  #number of cols equal?
  if len(vals1) != len(vals2): die("Number of columns different")
  
  #compare all cols with tolerance
  for j in range(0,len(vals1)):
    #if abs(float(vals1[j]) - float(vals2[j])) > options.eps: 
    diff=difference_relative(vals1[j], vals2[j], eps)
    if diff > eps: 
      die("mismatch in line %d col %d, %s != %s" %(i+1,j+1,vals1[j],vals2[j]))

