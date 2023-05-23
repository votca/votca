import unittest
import pytest
from ase.build import molecule
import os
try:
    from pyxtp import xtp
    from pyxtp.options import Options
except ImportError:
    raise unittest.SkipTest('pyxtp not available')


@pytest.fixture(scope="module")
def dftgwbse_xml():
    votcashare = os.environ.get("VOTCASHARE")
    if votcashare is None:
        msg = (
            "pyxtp: cannot find Votca installation, "
            "please set the VOTCASHARE environment variable"
        )
        raise RuntimeError(msg)
    return f"{votcashare}/xtp/xml/dftgwbse.xml"


class TestXTP:
    
    @pytest.fixture
    def CO(self):
        calculator = xtp(nthreads=2)
        atoms = molecule('CO', positions=([0,0,0],[1.4,0,0]), calculator=calculator)
        return atoms

    def test_get_and_read_potential_energy(self, CO):
        CO.get_potential_energy()
        new_atoms = xtp.read_atoms('CO.orb')
        assert(CO.calc.results['energy'] == new_atoms.calc.results['energy'])
         
    @pytest.mark.xfail()
    def test_failed_option(self, CO):
        CO.calc.options.dftpackage.basisset = 'invalid_basis_set'
        CO.get_potential_energy()
        
    def test_set(sefl, CO):
        CO.calc.set(basis='def2-svp', xc='lda', charge=1)
        assert CO.calc.options.dftpackage.basisset == 'def2-svp'
        assert CO.calc.options.dftpackage.functional == 'lda'
        assert CO.calc.options.dftpackage.charge == 1
        
    def test_set_default(self, CO):
        CO.calc.set(basis='sto-3g')
        CO.calc.set_default_options()
        assert CO.calc.options.dftpackage.basisset == ''
        
    def test_set_from_options(self, CO, dftgwbse_xml):
        opt = Options(dftgwbse_xml)
        opt.dftpackage.basisset = 'sto-3g'
        CO.calc.set_from_options(opt)
        assert CO.calc.options.dftpackage.basisset == 'sto-3g'
        
    @pytest.mark.xfail()
    def test_failed_get_total_energy(self, CO):
        CO.get_total_energy(name='invalid_energy_name')

    def test_get_total_energy(self, CO):
        CO.get_potential_energy()
        for name in CO.calc.implemented_properties:
            
            if name not in ["forces", "oscillator_strength"]:
                assert CO.calc.get_total_energy(name, level=1, dynamic = False) != None
                
            if name == 'oscillator_strength':
                energy, osc = CO.calc.get_oscillator_strengths(dynamic=False)
                assert energy.any() != None
                assert osc.any() != None
                
    def test_get_and_read_forces(self, CO):
        CO.calc.select_force(energy='singlets', level=1, dynamic=False)
        assert CO.get_forces().any() != None
             
        
        