from ase.utils import StringIO
from ase.io import read
from ase.io import xyz

# Made from NWChem interface


def read_geom_orcainp(filename):
    """Method to read geometry from an ORCA input file."""
    f = filename
    if isinstance(filename, str):
        f = open(filename)
    lines = f.readlines()
    if type(filename) == str:
        f.close()

    # Find geometry region of input file.
    stopline = 0
    for index, line in enumerate(lines):
        if line[1:].startswith('xyz '):
            startline = index + 1
            stopline = -1
        elif (line.startswith('end') and stopline == -1):
            stopline = index
        elif (line.startswith('*') and stopline == -1):
            stopline = index
    # Format and send to read_xyz.
    xyz_text = '%i\n' % (stopline - startline)
    xyz_text += ' geometry\n'
    for line in lines[startline:stopline]:
        xyz_text += line
    atoms = read(StringIO(xyz_text), format='xyz')
    atoms.set_cell((0., 0., 0.))  # no unit cell defined

    return atoms

