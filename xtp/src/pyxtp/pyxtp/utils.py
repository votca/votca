"""Utilities and constants."""
import scipy.constants
import numpy as np

H2EV = scipy.constants.physical_constants['Hartree energy in eV'][0]
BOHR2ANG = scipy.constants.physical_constants['Bohr radius'][0] * 1.e10
INVCM2EV = (scipy.constants.c *
            scipy.constants.physical_constants['Planck constant in eV/Hz'][0])*100.0
AFU2INVS = np.sqrt((scipy.constants.physical_constants['Hartree energy'][0]/np.power(
    scipy.constants.physical_constants['Bohr radius'][0], 2))/scipy.constants.physical_constants['atomic mass constant'][0])
AFU2INVCM = AFU2INVS/(2.0*np.pi*scipy.constants.c*100.)
AFU2EV = scipy.constants.physical_constants['reduced Planck constant in eV s'][0] * AFU2INVS
