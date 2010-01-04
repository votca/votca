#! /usr/bin/env python

#authors: James K and Bjorn Baumeier

import numpy as np
from cclib.parser import ccopen
import logging

# BEGIN changes BB 9-9-2009, read 'namejob' from command line call
# for calling the script from the command line, I add parsing of
# command line arguments:
# - the 'namejob' used for the dataset
# - the no. of the orbital on MolA 'orbA' and on MolB 'orbB' for which
#   the transfer integral is to be calculated (output)
#
# - calculated JAB_eff, instead of JAB
#  COMMAND LINE CALL:
#  ./ProJ-eff.py namejob orbA orbB



def usage():
    print ""
    print ("Script to determine the transfer integrals between two")
    print ("monomers by using a counterpoise-corrected projection")
    print ("on the molecular orbitals of a dimer.")
    print ("Call on command line:")
    print ""
    print ("  ./ProJ-eff.py namejob orbA orbB ")
    print ""
    print ("where:")
    print ("- namejob is the previously specified jobname")
    print ("- orbA    is the orbital of MolA, and")
    print ("- orbB    is the orbital of MolB")
    print ("for which the transfer integral is to be determined!")
    print ""
    print "Exiting!"
    exit(1)

import sys, getopt

# check if enough options are given on command line
if len(sys.argv)!=4:
    print("Error: Insufficient options given!")
    usage()

namejob=sys.argv[1]
orbA=sys.argv[2]
orbB=sys.argv[3]

# Changed arguments supplied to CalcJ to account for the specification
# or orbA and orbB on the command line
def CalcJ(namejob,orbA,orbB):
    molA_parser=ccopen("part1/"+namejob+"part1.log")
    molB_parser=ccopen("part2/"+namejob+"part2.log")
    molAB_parser=ccopen("dim/"+namejob+"dim.log")

    # Limit parsing info to ERROR messages
    molA_parser.logger.setLevel(logging.ERROR)
    molB_parser.logger.setLevel(logging.ERROR)
    molAB_parser.logger.setLevel(logging.ERROR)

    # Parse LOG files
    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()

    nbs=molAB.nbasis

    # Consistency checks
    if molA.nbasis!=molB.nbasis:
    print("Count of basis functions doesn't match. Failing.")
    return False

    for mole in molA,molB,molAB:
        if len(mole.atomcoords)!=1:
        print(mole+" calculation appears to be an optimisation! Failing.")

    #Take mocoeffs[0] - the alpha electrons to get matrix in correct order
    MolAB_Pro = np.transpose(np.dot(molAB.mocoeffs[0],molAB.aooverlaps))

    # Determine orbitals of monomers in the basis set of the dimer
    PsiA_DimBS = np.dot(molA.mocoeffs[0], MolAB_Pro)
    PsiB_DimBS = np.dot(molB.mocoeffs[0], MolAB_Pro)

    # Determine transfer integral JAB and site energies JAA, JBB
    JAB = np.dot(np.dot( PsiB_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiA_DimBS) )
    JAA = np.dot(np.dot( PsiA_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiA_DimBS) )
    JBB = np.dot(np.dot( PsiB_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiB_DimBS) )

    # Determine overlap analogous to JAB
    SAB = np.dot(PsiB_DimBS , np.transpose(PsiA_DimBS) )

    # Calculate JAB_eff according to Eq.10 in JACS 128, 9884 (2006)
    # !only for the desired orbitals!
    JAB_eff = (JAB[orbA,orbB] - 0.5*(JAA[orbA,orbA]+JBB[orbB,orbB])*SAB[orbA,orbB])/(1.0 - SAB[orbA,orbB]*SAB[orbA,orbB])

   

#    return JAB[orbA,orbB]
    return JAB_eff

# Finally, really call the CalcJ function and output the transfer integral
# for the orbitals specified by 'orbA' (on MolA) and 'orbB' (on MolB)
print CalcJ(namejob,orbA,orbB)

