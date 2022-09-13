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
from .options import Options
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
        self.options = Options(options)


    def update_options(self):
        """Merge user options with the defaults."""
        # parsing defaults
        votcashare = os.environ.get('VOTCASHARE')
        default_options_file = f'{votcashare}/xtp/xml/dftgwbse.xml'
        default_options = ET.parse(default_options_file)
        default_root = default_options.getroot()

        # prepare user options as ET
        user_options = ET.Element("options")
        user_options.append(create_xml_tree(ET.Element("dftgwbse"), self.options.to_dict()))
        ET.ElementTree(user_options).write("test.xml")
        ET.ElementTree(user_options).write('dftgwbse.xml')
  

    def run(self, nThreads: int, path_examples: Path):
        """Just runs xtp_tools with command line call."""
        # update and write the options
        self.update_options()

        # write the XYZfile
        xyzname = self.mol.name
        xyzfile = xyzname + ".xyz"
        self.mol.write_xyz_file(xyzfile)

        """ Runs VOTCA and moves results a job folder, if requested """
        if not Path(self.jobdir).exists():
            os.makedirs(self.jobdir)

        path_dftgwbse = (path_examples / "dftgwbse.xml").absolute().as_posix()
        
        # Call the tool and capture the standard output
        output = capture_standard_output(
            xtp_binds.call_tool, "dftgwbse", nThreads, path_dftgwbse)

        with open("example.out", "w") as handler:
            handler.write(output)



        # votcacmd = f"xtp_tools -e dftgwbse -c job_name={xyzname} -o dftgwbse.xml -t {self.threads} > {self.jobdir}{self.jobname}.log"

        # # DYLD_LIBRARY_PATH is removed from env under MacOS
        # if platform.system() == "Darwin":
        #     DYLD=os.environ.get('DYLD_LIBRARY_PATH')
        #     votcacmd = f"export DYLD_LIBRARY_PATH={DYLD};" + votcacmd
        
        # subprocess.run(votcacmd, shell=True, stderr=subprocess.STDOUT)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = f'{self.jobdir}{self.jobname}.orb'
            os.replace(f'{xyzname}.orb', self.orbfile)
        else:
            self.orbfile = f'{xyzname}.orb'

        # Reads energies from an existing HDF5 orb file
        self.mol.read_orb(self.orbfile)

