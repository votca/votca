from pyxtp import Molecule
from pathlib import Path
import os

PATH_TEST = Path(os.getcwd()) / "test/files"

def test_molecule_io(tmp_path: Path):
    """Check molecule methods."""
    # check reading method
    mol = Molecule.from_xyz_file(PATH_TEST / "ethylene.xyz")

    assert len(mol.elements) == 6

    # Check writing method
    file_name = tmp_path / "test.xyz"
    Molecule.to_xyz_file(file_name, mol)
    assert file_name.exists()
