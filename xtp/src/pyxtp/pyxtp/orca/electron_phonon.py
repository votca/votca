"""Electron-phonon-coupling module."""

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import periodictable
import scipy.linalg as la

from ..utils import AFU2EV, H2EV


class Electronphonon:
    
    def get_mass(self, elements):
        mass = []
        for atom in elements:
            mass.append(periodictable.formula(atom).mass)
        mass = np.array(mass, dtype=float)
        return mass

    def calculate_electron_phonon_couplings(
            self, elements: np.ndarray, gs_hessian: np.ndarray, gradient: np.ndarray,
            plot: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """Calculates the mode-resolved electron-phonon coupling from 
           the normal mode projected excited state gradient:

                V_nu^ep = 1/(M * w_nu^2) * | e_v * gradient(E) |^2

            where:
                - w_nu: frequency of the nu-th vibrational mode in GS
                - e_nu: eigenvector of the nu-th vibrational mode in GS

            Parameters
            ----------
            elements  np.array with N chemical elements 
            gs_hessian  np.array of GS Hessian (3N,3N) in (a.f.u)^2 = Hartree/Bohr^2/amu
            gradient  flattened excited state gradient (3N) in Hartree/Bohr

            Returns
            -------
            freq vibrational frequencies in a.f.u
            ep_couplings electron-phonon coupling per mode in Hartree

        """

        # we don't want to change the input gradient
        es_gradient = np.copy(gradient)

        # get vector with masses of the atoms
        mass = np.repeat(self.get_mass(elements), 3)

        # get vibrational frequencies and eigenmodes
        freq, modes = self.calculate_vibrational_modes(elements, gs_hessian)

        # flatten the excited state gradient, if necessary
        es_gradient = es_gradient.flatten()

        # get the normal-mode projected gradient (already divided by sqrt(M))
        es_gradient /= np.sqrt(mass)
        nm_gradient = modes.T.dot(es_gradient)

        # ignore any contribution from translation/rotation and imaginary frequencies
        # and frequencies lower than 0.002 a.f.u ~ 10cm^-1 
        nm_gradient[np.where(freq < 2e-3)] = 0.

        # electron-phonon couplings in Hartree
        ep_couplings = (nm_gradient / freq) ** 2

        if plot:
            plt.bar(AFU2EV * freq, H2EV * ep_couplings, width=0.005, alpha=0.6)
            plt.plot(AFU2EV * freq, H2EV * np.cumsum(ep_couplings), lw=2., alpha=0.8)
            plt.ylabel('cum. electron-phonon coupling (eV)')
            plt.xlabel('hw (eV)')
            plt.show()

        return (freq, ep_couplings)

    def calculate_vibrational_modes(self, elements, hessian_in) -> Tuple[np.ndarray, np.ndarray]:
        """Calculates the vibrational frequencies and eigenmodes from a Hessian matrix and a list of elements"""

        # we don't want to change the input hessian
        hessian = np.copy(hessian_in)

        # get vector with masses of the atoms
        mass = np.repeat(self.get_mass(elements), 3)

        # mass weighting the Hessian
        hessian /= np.sqrt(mass)
        hessian /= np.sqrt(mass.reshape((len(mass), 1)))

        # get eigenmodes (squared frequency in (a.f.u)^2 = Hartree/Bohr^2/amu)
        freq_sq, modes = la.eig(hessian)

        # convert to frequencies in a.f.u
        cfreq = np.array(freq_sq, dtype=np.complex128)
        cfreq = np.sqrt(cfreq)

        # check for imaginary frequencies and store them as negative real ones
        freq = np.where(np.isreal(cfreq), np.real(cfreq), -np.abs(np.imag(cfreq)))

        # sort frequencies and eigenvectors
        isort = freq.argsort()
        modes = np.real(modes[:, isort])
        freq = np.sort(freq)

        return (freq, modes)
