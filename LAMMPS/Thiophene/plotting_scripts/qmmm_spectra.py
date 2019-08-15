import h5py
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--orbfile",required=True, help="checkpoint HDF5 file",)
parser.add_argument("-o", "--outputfile",default="spectrum.dat", help="output file for spectra")
args = parser.parse_args()

def getEnergies(orb):
    a=orb['region_0']['orbitals']['BSE_singlet']['eigenvalues'].value
    a.flatten()
    return a.flatten()

def trans_sort(index):
    return int(index[3:])

def getOscillators(orb):
    energies=getEnergies(orb)
    transdip=[]
    for k in sorted(orb['region_0']['orbitals']['transition_dipoles'].keys(),key=trans_sort):
        print(k)
        transdip.append(np.array(orb['region_0']['orbitals']['transition_dipoles'][k].value)) 
    d2=[]
    for b in transdip:
        d2.append(np.sum(b**2))    
    d2=np.array(d2)
    oscs=2/3.0*energies*d2
    return oscs

orb=h5py.File(args.orbfile,'r')
e=getEnergies(orb)*27.2114
osc=getOscillators(orb)
output=np.array([e,osc])
np.savetxt(args.outputfile,output.T,header="Omega[eV],f")
print("Wrote spectrum to {}.".format(args.outputfile))
