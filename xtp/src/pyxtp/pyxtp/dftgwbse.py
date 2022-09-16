"""DFTGWSE wrapper."""
import os
import platform
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional

from .capture_standard_output import capture_standard_output
from pyxtp import xtp_binds

from .molecule import Molecule
from .options import Options, XTPOptions
from .xml_editor import create_xml_tree

__all__ = ["DFTGWBSE"]


class DFTGWBSE:

    def __init__(self, mol: Molecule, threads: int = 1, jobname: str = 'dftgwbse',
                 options: Optional[Dict[str, Any]] = {}, jobdir: str = './'):
        self.mol = mol
        self.threads = threads
        self.jobname = jobname
        self.jobdir = jobdir
        self.orbfile = ''
        self.options = XTPOptions()

  
    def run(self, nThreads: int = 1):
        """Just runs xtp_tools with command line call."""

        # write the XYZfile
        xyzname = self.mol.name
        xyzfile = xyzname + ".xyz"
        self.mol.write_xyz_file(xyzfile)

        # update and write the options
        self.options.data.job_name = xyzname
        # self.options.write_xml()
        # self.update_options()

        """ Runs VOTCA and moves results a job folder, if requested """
        if not Path(self.jobdir).exists():
            os.makedirs(self.jobdir)

        # path_dftgwbse = (path_examples / "dftgwbse.xml").absolute().as_posix()
        # path_dftgwbse = (Path("./") / "dftgwbse.xml").absolute().as_posix()
        #path_dftgwbse = (Path("files_examples") / "dftgwbse.xml").absolute().as_posix()
        path_dftgwbse = (Path("dftgwbse.xml")).absolute().as_posix()
        
        # Call the tool and capture the standard output
        output = capture_standard_output(
            xtp_binds.call_tool, "dftgwbse", nThreads, path_dftgwbse)

        with open("example.out", "w") as handler:
            handler.write(output)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = f'{self.jobdir}{self.jobname}.orb'
            os.replace(f'{xyzname}.orb', self.orbfile)
        else:
            self.orbfile = f'{xyzname}.orb'

        # Reads energies from an existing HDF5 orb file
        self.mol.read_orb(self.orbfile)

