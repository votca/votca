#python script to perform singular value decompositon of the cross correlation matrix
#copy the file in the folder of the first imc iteration and run: python svd.py

#!/bin/env python2

import numpy as np
import numpy.linalg as la
np.set_printoptions(threshold=np.nan)
A=np.loadtxt('group_1.gmc');
b=np.loadtxt('group_1.imc');


#singular value decomposition
[U,s,V] =la.svd(A);
#write singular values
np.savetxt('singular_values.dat',s)
