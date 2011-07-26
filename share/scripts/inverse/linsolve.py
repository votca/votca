#!/usr/bin/env python

import numpy as np
A = np.loadtxt('$name.gmc');
b = np.loadtxt('$name.imc');

x = np.empty([len(b),2])
x[:,0]=b[:,0];
x[:,1] = -np.linalg.solve(A,b[:,1]);

np.savetxt('$name_out', x)
