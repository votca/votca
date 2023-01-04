"""DFTGWSE wrapper."""
import os
from pathlib import Path
from typing import Any, Dict, Optional

from .capture_standard_output import capture_standard_output
from pyxtp import xtp_binds

from .molecule import Molecule
from .options import XTPOptions

__all__ = ["DFTGWBSE"]


class DFTGWBSE:

    def __init__(self, mol: Molecule):
        self.mol = mol
        self.orbfile = ''
        self.options = XTPOptions()
        self.jobdir = "./"

  
    def run(self, nThreads: int = 1):
        """Just runs xtp_tools with command line call."""

        # write the XYZfile
        xyzname = self.mol.name
        xyzfile = xyzname + ".xyz"
        self.mol.write_xyz_file(xyzfile)

        # update and write the options
        self.options.job_name = xyzname
        input_filename='dftgwbse.xml'
        self.options._write_xml(input_filename)


        """ Runs VOTCA and moves results a job folder, if requested """
        if not Path(self.jobdir).exists():
            os.makedirs(self.jobdir)

        # current dftgwbse file
        path_dftgwbse = (Path(input_filename)).absolute().as_posix()
        
        # Call the tool and capture the standard output
        output = capture_standard_output(
            xtp_binds.call_tool, "dftgwbse", nThreads, path_dftgwbse)

        with open(xyzname + ".out", "w") as handler:
            handler.write(output)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = f'{self.jobdir}{self.options.job_name}.orb'
            os.replace(f'{xyzname}.orb', self.orbfile)
        else:
            self.orbfile = f'{xyzname}.orb'

        # Reads energies from an existing HDF5 orb file
        self.mol.read_orb(self.orbfile)

