"""Visualization module."""
from .molecule import Molecule
import numpy as np
import matplotlib.pyplot as plt
from .utils import H2EV


class Visualization:
    def __init__(self, mol: Molecule, save_figure: bool = False):
        self.mol = mol
        self.save_figure = save_figure

    def plot_qp_corrections(self):
        qp_corrections = self.mol.get_qp_corrections()
        qpmin = self.mol.qpmin
        qpmax = self.mol.qpmax + 1
        corrected_ks = self.mol.ks_energies[qpmin:qpmax]
        plt.plot(corrected_ks, qp_corrections, 'ro')
        plt.xlabel('KS Energy (eV)')
        plt.ylabel('QP correction (eV)')
        if self.save_figure:
            plt.savefig("qpcorrections.png")
        else:
            plt.show()

    def plot_absorption_gaussian(
            self, dynamic: bool = False, min: float = 0.0, max: float = 10.0, points: int = 1000,
            sigma: float = 0.2):

        energy, osc = self.mol.get_oscillator_strengths(dynamic)

        # convert energies from Hartree to eV
        energy *= H2EV  # to eV
        # make a stick plot with oscillator strength
        plt.stem(energy, osc, basefmt=" ")
        # apply unormalized Gaussian lineshape
        e = np.linspace(min, max, points)
        spectrum = 0

        for i in range(len(energy)):
            spectrum += osc[i] * self.gaussian(e, energy[i], sigma)

        plt.plot(e, spectrum, 'k', linewidth=2)
        plt.ylim(bottom=0)
        plt.title(f'Gaussian lineshape with sigma = {sigma}eV')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Absorption (arb. units)')
        if self.save_figure:
            plt.savefig("absorption_gaussian.png")
        else:
            plt.show()


    def gaussian(self, x, mu, sig):
        """Non-normalized Gaussian distribution."""
        return np.exp(-0.5 * ((x - mu) / sig) ** 2)

    def plot_dos_gaussian(self, ks=True, qp=True, min: float = -25.0, max: float = 10.0, points: int = 1000,
                          sigma: float = 0.2):

        plt.xlabel('Energy (eV)')
        plt.ylabel('DOS (arb. units)')
        e = np.linspace(min, max, points)
        dos_ks = 0
        if ks:
            ks_energies = self.mol.ks_energies * H2EV  # to eV
            qpmin = self.mol.qpmin
            qpmax = self.mol.qpmax + 1
            for ks_energy in ks_energies[qpmin:qpmax]:
                dos_ks += self.gaussian(e, ks_energy, sigma)
            plt.plot(e, dos_ks, 'k', linewidth=2, label='KS')
        dos_qp = 0
        if qp:
            qp_energies = self.mol.qp_energies * H2EV 
            for qp_energy in qp_energies:
                dos_qp += self.gaussian(e, qp_energy, sigma)
            plt.plot(e, dos_qp, 'r', linewidth=2, label='QP')
        plt.legend()
        if self.save_figure:
            plt.savefig("dos_gaussian.png")
        else:
            plt.show()
