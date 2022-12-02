import os
import numpy as np
import h5py

from ase.units import Hartree, Bohr
from ase.io.votca import write_votca
from ase.calculators.calculator import FileIOCalculator, Parameters, ReadError


class votca(FileIOCalculator):
    implemented_properties = ['energy', 'forces', 'singlets',
        'triplets', 'qp', 'ks', 'qp_pert', 'transition_dipoles']

    command = 'xtp_tools -e dftgwbse -o dftgwbse.xml -t 4 >  dftgwbse.log'

    default_parameters = dict(
        charge=0, mult=1,
        task='gradient',
        orcasimpleinput='tightscf PBE def2-SVP',
        orcablocks='%scf maxiter 200 end')

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='orca', atoms=None, **kwargs):
        """ ASE interface to VOTCA-XTP
        Only supports energies for now.
        """
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.pcpot = None

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        p = self.parameters
        p.write(self.label + '.ase')
        p['label'] = self.label
        # if self.pcpot:  # also write point charge file and add things to input
        #    p['pcpot'] = self.pcpot

        write_votca(atoms, **p)

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile(self.label + '.out'):
            raise ReadError

        with open(self.label + '.inp') as f:
            for line in f:
                if line.startswith('geometry'):
                    break
            symbols = []
            positions = []
            for line in f:
                if line.startswith('end'):
                    break
                words = line.split()
                symbols.append(words[0])
                positions.append([float(word) for word in words[1:]])

        self.parameters = Parameters.read(self.label + '.ase')
        self.read_results(self)

    def read_results(self):
        self.read_energy()
        if self.parameters.task.find('gradient') > -1:
            self.read_forces()

    def tdipoles_sorter(orb):
      groupedData = []
      permutation = []
      for ind in orb['transition_dipoles'].keys():
            permutation.append(int(ind[3:]))  # 3: skips over the 'ind' bit
            groupedData.append(orb['transition_dipoles'][ind][:].transpose()[0])
      groupedData = np.asarray(groupedData)
      return(groupedData[np.argsort(permutation)])

    def read_energy(self):
        """Read Energy from VOTCA-XTP log file."""
        orbFile = h5py.File('system.orb', 'r')
        orb = orbFile['QMdata']
        self.results['energy'] = orb.attrs['qm_energy']
        self.results['singlets'] = np.array(
            orb['BSE_singlet']['eigenvalues'][()]).transpose()[0]
        self.results['triplets'] = np.array(
            orb['BSE_triplet']['eigenvalues'][()]).transpose()[0]
        self.results['ks'] = np.array(
            orb['mos']['eigenvalues'][()]).transpose()[0]
        self.results['qp'] = np.array(
            orb['QPdiag']['eigenvalues'][()]).transpose()[0]
        groupedData = []
        permutation = []
        for ind in orb['transition_dipoles'].keys():
            permutation.append(int(ind[3:])) # 3: skips over the 'ind' bit
            groupedData.append(orb['transition_dipoles'][ind][:].transpose()[0])
        groupedData = np.asarray(groupedData)
        self.results['transition_dipoles'] = groupedData[np.argsort(permutation)]
        self.results['qp_pert'] =  np.array(orb['QPpert_energies'][()]).transpose()[0]
    def read_forces(self):
        """Read Forces from VOTCA logfile."""
        fil = open('dftgwbse.log', 'r', encoding='utf-8')
        lines = fil.readlines()
        fil.close()
        getgrad = "no"
        for i, line in enumerate(lines):
            if line.find('ENGRAD') >= 0:
                getgrad = "yes"
                gradients = []
                tempgrad = []
                continue
            if getgrad == "yes" and "Saving" not in line:
                if "Total" not in line:
                    grad = line.split()
                    tempgrad.append(float(grad[3]))
                    tempgrad.append(float(grad[4]))
                    tempgrad.append(float(grad[5]))
                    gradients.append(tempgrad)
                    tempgrad = []
                else:
                    energy = float(line.split()[-2])
            if 'Saving' in line:
                getgrad = "no"
        self.results['forces'] = -np.array(gradients) * Hartree / Bohr

    def embed(self, mmcharges=None, **parameters):
        """Embed atoms in point-charges (mmcharges)
        """
        self.pcpot = PointChargePotential(mmcharges, label=self.label)
        return self.pcpot


