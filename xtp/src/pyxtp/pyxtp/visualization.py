"""Visualization module."""
from ase import Atoms
import numpy as np
import matplotlib.pyplot as plt
from .utils import H2EV


class Visualization:
    def __init__(self, atoms: Atoms, save_figure: bool = False):
        self.atoms = atoms
        self.save_figure = save_figure

    def plot_qp_corrections(self):
        qp_corrections = self.atoms._calc.get_qp_corrections()
        qpmin = self.atoms._calc.qpmin
        qpmax = self.atoms._calc.qpmax + 1
        corrected_ks = self.atoms._calc.ks_energies[qpmin:qpmax]
        plt.plot(corrected_ks, qp_corrections, 'ro')
        plt.xlabel('KS Energy (eV)')
        plt.ylabel('QP correction (eV)')
        if self.save_figure:
            plt.savefig("qpcorrections.png")
        else:
            plt.show()

    def plot_absorption_gaussian(
            self, dynamic: bool = False, energy_min: float = 0.0, energy_max: float = 10.0, points: int = 1000,
            sigma: float = 0.2):

        energy, osc = self.atoms._calc.get_oscillator_strengths(dynamic)

        # convert energies from Hartree to eV
        energy *= H2EV  # to eV
        # make a stick plot with oscillator strength
        plt.stem(energy, osc, basefmt=" ")
        # apply unormalized Gaussian lineshape
        e = np.linspace(energy_min, energy_max, points)
        spectrum = 0

        for i, eval in enumerate(energy):
            spectrum += osc[i] * self.gaussian(e, eval, sigma)

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

    def plot_dos_gaussian(self, ks=True, qp=True, energy_min: float = -25.0, energy_max: float = 10.0, points: int = 1000,
                          sigma: float = 0.2):

        plt.xlabel('Energy (eV)')
        plt.ylabel('DOS (arb. units)')
        e = np.linspace(energy_min, energy_max, points)
        dos_ks = 0
        if ks:
            ks_energies = self.atoms._calc.ks_energies * H2EV  # to eV
            qpmin = self.atoms._calc.qpmin
            qpmax = self.atoms._calc.qpmax + 1
            for ks_energy in ks_energies[qpmin:qpmax]:
                dos_ks += self.gaussian(e, ks_energy, sigma)
            plt.plot(e, dos_ks, 'k', linewidth=2, label='KS')
        dos_qp = 0
        if qp:
            qp_energies = self.atoms._calc.qp_energies * H2EV
            for qp_energy in qp_energies:
                dos_qp += self.gaussian(e, qp_energy, sigma)
            plt.plot(e, dos_qp, 'r', linewidth=2, label='QP')
        plt.legend()
        if self.save_figure:
            plt.savefig("dos_gaussian.png")
        else:
            plt.show()
