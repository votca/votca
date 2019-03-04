#!/ usr / bin / env python3
#
#Copyright 2009 - 2019 The VOTCA Development Team(http:  // www.votca.org)
#
#Licensed under the Apache License, Version 2.0(the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#http:  // www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#

from optparse import OptionParser import numpy as np import numpy.linalg as la

    usage = "Usage: %prog [options] group output" parser =
    OptionParser(usage = usage) parser.add_option(
        "--reg", dest = "reg", metavar = "REG", help = "regularization factor",
        default = 0)(options, args) =
        parser.parse_args()

            if len(args) != 2 : exit("two statefile required as parameters")

                                    A = np.loadtxt(args[0] + '.gmc');
b                                     = np.loadtxt(args[0] + '.imc');
x                                     = np.empty([ len(b), 2 ]) n,
    m = A.shape I  = np.identity(
        m) x[:, 0] = b
            [:, 0] x
            [:, 1] = -np.dot(
                       np.dot(la.inv(np.dot(A.T, A) + float(options.reg) * I),
                              A.T),
                       b[:, 1]);

np.savetxt(args[1], x)
