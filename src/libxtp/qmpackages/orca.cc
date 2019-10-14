/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "orca.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <iomanip>
#include <stdio.h>
#include <votca/tools/elements.h>
#include <votca/xtp/ecpaobasis.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
using namespace std;

void Orca::Initialize(tools::Property& options) {

  // good luck

  // Orca file names
  std::string fileName = "system";

  _input_file_name = fileName + ".inp";
  _log_file_name = fileName + ".log";
  _shell_file_name = fileName + ".sh";
  _mo_file_name = fileName + ".gbw";

  ParseCommonOptions(options);

  // check if the optimize keyword is present, if yes, read updated coords
  std::string::size_type iop_pos =
      _options.find(" Opt"); /*optimization word in orca*/
  if (iop_pos != std::string::npos) {
    _is_optimization = true;
  }

  if (_write_guess) {
    iop_pos = _options.find("Guess MORead");
    if (iop_pos != std::string::npos) {
      _options = _options + "\n Guess MORead ";
    }
  }
}

/* Custom basis sets are written on a per-element basis to
 * the system.bas/aux file(s), which are then included in the
 * Orca input file using GTOName = "system.bas/aux"
 */
void Orca::WriteBasisset(const QMMolecule& qmatoms, std::string& bs_name,
                         std::string& el_file_name) {

  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  tools::Elements elementInfo;
  BasisSet bs;
  bs.Load(bs_name);
  XTP_LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << bs_name << flush;
  ofstream el_file;

  el_file.open(el_file_name);
  el_file << "$DATA" << endl;

  for (const std::string& element_name : UniqueElements) {
    const Element& element = bs.getElement(element_name);
    el_file << elementInfo.getEleFull(element_name) << endl;
    for (const Shell& shell : element) {
      for (const char& subtype : shell.getType()) {
        el_file << subtype << " " << shell.getSize() << endl;
        int sh_idx = 0;
        for (const GaussianPrimitive& gaussian : shell) {
          sh_idx++;
          el_file << " " << sh_idx << " " << indent(gaussian.decay());
          el_file << " "
                  << indent(gaussian.Contractions()[FindLmax(
                         std::string(1, subtype))]);
          el_file << endl;
        }
      }
    }
  }
  el_file << "STOP\n";
  el_file.close();

  return;
}

/* Coordinates are written in standard Element,x,y,z format to the
 * input file.
 */
void Orca::WriteCoordinates(std::ofstream& inp_file,
                            const QMMolecule& qmatoms) {

  for (const QMAtom& atom : qmatoms) {
    Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    inp_file << setw(3) << atom.getElement() << setw(12)
             << setiosflags(ios::fixed) << setprecision(5) << pos.x()
             << setw(12) << setiosflags(ios::fixed) << setprecision(5)
             << pos.y() << setw(12) << setiosflags(ios::fixed)
             << setprecision(5) << pos.z() << endl;
  }
  inp_file << "* \n" << endl;
  return;
}

/* If custom ECPs are used, they need to be specified in the input file
 * in a section following the basis set includes.
 */
void Orca::WriteECP(std::ofstream& inp_file, const QMMolecule& qmatoms) {

  inp_file << endl;
  std::vector<std::string> UniqueElements = qmatoms.FindUniqueElements();

  ECPBasisSet ecp;
  ecp.Load(_ecp_name);

  XTP_LOG(logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

  for (const std::string& element_name : UniqueElements) {
    try {
      ecp.getElement(element_name);
    } catch (std::runtime_error& error) {
      XTP_LOG(logDEBUG, *_pLog)
          << "No pseudopotential for " << element_name << " available" << flush;
      continue;
    }
    const ECPElement& element = ecp.getElement(element_name);

    inp_file << "\n"
             << "NewECP"
             << " " << element_name << endl;
    inp_file << "N_core"
             << " " << element.getNcore() << endl;
    inp_file << "lmax"
             << " " << getLName(element.getLmax()) << endl;
    // For Orca the order doesn't matter but let's write it in ascending order
    // write remaining shells in ascending order s,p,d...
    for (int i = 0; i <= element.getLmax(); i++) {
      for (const ECPShell& shell : element) {
        if (shell.getL() == i) {
          // shell type, number primitives, scale factor
          inp_file << shell.getType() << " " << shell.getSize() << endl;
          int sh_idx = 0;
          for (const ECPGaussianPrimitive& gaussian : shell) {
            sh_idx++;
            inp_file << sh_idx << " " << gaussian._decay << " "
                     << gaussian._contraction << " " << gaussian._power << endl;
          }
        }
      }
    }
    inp_file << "end\n "
             << "\n"
             << endl;
  }
  return;
}

void Orca::WriteChargeOption() {
  std::string::size_type iop_pos = _options.find("pointcharges");
  if (iop_pos == std::string::npos) {
    _options = _options + "\n %pointcharges \"background.crg\"";
  }
}

/* For QM/MM the molecules in the MM environment are represented by
 * their atomic partial charge distributions. ORCA expects them in
 * q,x,y,z format in a separate file "background.crg"
 */
void Orca::WriteBackgroundCharges() {

  std::ofstream crg_file;
  std::string _crg_file_name_full = _run_dir + "/background.crg";
  crg_file.open(_crg_file_name_full);
  int total_background = 0;

  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    if (site->getCharge() != 0.0) total_background++;
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    total_background += split_multipoles.size();
  }  // counting only

  crg_file << total_background << endl;
  boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
  // now write
  for (const std::unique_ptr<StaticSite>& site : _externalsites) {
    Eigen::Vector3d pos = site->getPos() * tools::conv::bohr2ang;
    string sitestring =
        boost::str(fmt % site->getCharge() % pos.x() % pos.y() % pos.z());
    if (site->getCharge() != 0.0) crg_file << sitestring << endl;
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    for (const auto& mpoles : split_multipoles) {
      Eigen::Vector3d pos = mpoles._pos * tools::conv::bohr2ang;
      string multipole =
          boost::str(fmt % mpoles._q % pos.x() % pos.y() % pos.z());
      crg_file << multipole << endl;
    }
  }

  return;
}

/**
 * Prepares the *.inp file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied, Not
 * implemented yet
 */
bool Orca::WriteInputFile(const Orbitals& orbitals) {

  std::vector<std::string> results;
  std::string temp_suffix = "/id";
  std::string scratch_dir_backup = _scratch_dir;
  std::ofstream inp_file;
  std::string inp_file_name_full = _run_dir + "/" + _input_file_name;
  inp_file.open(inp_file_name_full);
  // header
  inp_file << "* xyz  " << _charge << " " << _spin << endl;
  int threads = OPENMP::getMaxThreads();
  const QMMolecule& qmatoms = orbitals.QMAtoms();
  // put coordinates
  WriteCoordinates(inp_file, qmatoms);
  // add parallelization info
  inp_file << "%pal\n "
           << "nprocs " << threads << "\nend"
           << "\n"
           << endl;
  // basis set info
  if (_write_basis_set) {
    std::string el_file_name = _run_dir + "/" + "system.bas";
    WriteBasisset(qmatoms, _basisset_name, el_file_name);
    inp_file << "%basis\n " << endl;
    inp_file << "GTOName"
             << " "
             << "="
             << "\"system.bas\";" << endl;
    if (_write_auxbasis_set) {
      std::string aux_file_name = _run_dir + "/" + "system.aux";
      WriteBasisset(qmatoms, _auxbasisset_name, aux_file_name);
      inp_file << "GTOAuxName"
               << " "
               << "="
               << "\"system.aux\";" << endl;
    }
  }  // write_basis set

  // ECPs
  if (_write_pseudopotentials) {
    WriteECP(inp_file, qmatoms);
  }
  inp_file << "end\n "
           << "\n"
           << endl;  // This end is for the basis set block
  if (_write_charges) {
    WriteBackgroundCharges();
  }

  inp_file << _options << "\n";
  inp_file << endl;
  inp_file.close();
  // and now generate a shell script to run both jobs, if neccessary

  XTP_LOG(logDEBUG, *_pLog)
      << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;
  _scratch_dir = scratch_dir_backup + temp_suffix;
  WriteShellScript();
  _scratch_dir = scratch_dir_backup;
  return true;
}

bool Orca::WriteShellScript() {
  ofstream shell_file;
  std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;
  shell_file.open(shell_file_name_full);
  shell_file << "#!/bin/bash" << endl;
  shell_file << "mkdir -p " << _scratch_dir << endl;

  if (_write_guess) {
    if (!(boost::filesystem::exists(_run_dir + "/molA.gbw") &&
          boost::filesystem::exists(_run_dir + "/molB.gbw"))) {
      throw runtime_error(
          "Using guess relies on a molA.gbw and a molB.gbw file being in the "
          "directory.");
    }
    shell_file << _executable
               << "_mergefrag molA.gbw molB.gbw dimer.gbw > merge.log" << endl;
  }
  shell_file << _executable << " " << _input_file_name << " > "
             << _log_file_name << endl;  //" 2> run.error" << endl;
  shell_file.close();
  return true;
}

/**
 * Runs the Orca job.
 */
bool Orca::Run() {

  XTP_LOG(logDEBUG, *_pLog) << "Running Orca job" << flush;

  if (std::system(NULL)) {

    std::string command = "cd " + _run_dir + "; sh " + _shell_file_name;
    int check = std::system(command.c_str());
    if (check == -1) {
      XTP_LOG(logERROR, *_pLog)
          << _input_file_name << " failed to start" << flush;
      return false;
    }
    if (CheckLogFile()) {
      XTP_LOG(logDEBUG, *_pLog) << "Finished Orca job" << flush;
      return true;
    } else {
      XTP_LOG(logDEBUG, *_pLog) << "Orca job failed" << flush;
    }
  } else {
    XTP_LOG(logERROR, *_pLog)
        << _input_file_name << " failed to start" << flush;
    return false;
  }

  return true;
}

/**
 * Cleans up after the Orca job
 */
void Orca::CleanUp() {

  if (_write_guess) {
    remove((_run_dir + "/" + "molA.gbw").c_str());
    remove((_run_dir + "/" + "molB.gbw").c_str());
    remove((_run_dir + "/" + "dimer.gbw").c_str());
  }
  // cleaning up the generated files
  if (_cleanup.size() != 0) {
    tools::Tokenizer tok_cleanup(_cleanup, ",");
    std::vector<std::string> cleanup_info;
    tok_cleanup.ToVector(cleanup_info);
    for (const std::string& substring : cleanup_info) {
      if (substring == "inp") {
        std::string file_name = _run_dir + "/" + _input_file_name;
        remove(file_name.c_str());
      }

      if (substring == "bas") {
        std::string file_name = _run_dir + "/system.bas";
        remove(file_name.c_str());
      }

      if (substring == "log") {
        std::string file_name = _run_dir + "/" + _log_file_name;
        remove(file_name.c_str());
      }

      if (substring == "gbw") {
        std::string file_name = _run_dir + "/" + _mo_file_name;
        remove(file_name.c_str());
      }

      if (substring == "ges") {
        std::string file_name = _run_dir + "/system.ges";
        remove(file_name.c_str());
      }
      if (substring == "prop") {
        std::string file_name = _run_dir + "/system.prop";
        remove(file_name.c_str());
      }
    }
  }
  return;
}

StaticSegment Orca::GetCharges() const {

  StaticSegment result("charges", 0);

  XTP_LOG(logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
  std::string log_file_name_full = _run_dir + "/" + _log_file_name;
  std::string line;

  std::ifstream input_file(log_file_name_full);
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);
    GetCoordinates(result, line, input_file);

    std::string::size_type charge_pos = line.find("CHELPG Charges");

    if (charge_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting charges" << flush;
      getline(input_file, line);
      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      int nfields = row.size();
      bool hasAtoms = result.size() > 0;
      while (nfields == 4) {
        int atom_id = boost::lexical_cast<int>(row.at(0));
        std::string atom_type = row.at(1);
        double atom_charge = boost::lexical_cast<double>(row.at(3));
        row = GetLineAndSplit(input_file, "\t ");
        nfields = row.size();
        if (hasAtoms) {
          StaticSite& temp = result.at(atom_id);
          if (temp.getElement() != atom_type) {
            throw std::runtime_error(
                "Getting charges failed. Mismatch in elemts:" +
                temp.getElement() + " vs " + atom_type);
          }
          temp.setCharge(atom_charge);
        } else {
          StaticSite temp =
              StaticSite(atom_id, atom_type, Eigen::Vector3d::Zero());
          temp.setCharge(atom_charge);
          result.push_back(temp);
        }
      }
    }
  }
  return result;
}

Eigen::Matrix3d Orca::GetPolarizability() const {
  std::string line;
  ifstream input_file((_run_dir + "/" + _log_file_name));
  bool has_pol = false;

  Eigen::Matrix3d pol = Eigen::Matrix3d::Zero();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    std::string::size_type pol_pos = line.find("THE POLARIZABILITY TENSOR");
    if (pol_pos != std::string::npos) {
      XTP_LOG(logDEBUG, *_pLog) << "Getting polarizability" << flush;
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);

      if (line.find("The raw cartesian tensor (atomic units)") ==
          std::string::npos) {
        throw std::runtime_error(
            "Could not find cartesian polarisation tensor");
      }

      for (int i = 0; i < 3; i++) {
        getline(input_file, line);
        tools::Tokenizer tok2(line, " ");
        std::vector<std::string> values = tok2.ToVector();
        if (values.size() != 3) {
          throw std::runtime_error("Polarisation line " + line +
                                   " cannot be parsed");
        }
        Eigen::Vector3d row;
        row << std::stod(values[0]), std::stod(values[1]), std::stod(values[2]);
        pol.row(i) = row;
      }

      has_pol = true;
    }
  }
  if (!has_pol) {
    throw std::runtime_error("Could not find polarisation in logfile");
  }
  return pol;
}

bool Orca::ParseLogFile(Orbitals& orbitals) {
  bool found_success = false;
  orbitals.setQMpackage(getPackageName());
  orbitals.setDFTbasisName(_basisset_name);
  if (_write_pseudopotentials) {
    orbitals.setECPName(_ecp_name);
  }

  XTP_LOG(logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
  std::string log_file_name_full = _run_dir + "/" + _log_file_name;
  // check if LOG file is complete
  if (!CheckLogFile()) return false;
  std::map<int, double> energies;
  std::map<int, double> occupancy;

  std::string line;
  int levels = 0;
  int number_of_electrons = 0;
  std::vector<std::string> results;

  std::ifstream input_file(log_file_name_full);

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog)
        << "File " << log_file_name_full << " not found " << flush;
    return false;
  } else {
    XTP_LOG(logDEBUG, *_pLog)
        << "Reading Coordinates and occupationnumbers and energies from "
        << log_file_name_full << flush;
  }
  // Coordinates of the final configuration depending on whether it is an
  // optimization or not

  QMMolecule& mol = orbitals.QMAtoms();
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);

    GetCoordinates(mol, line, input_file);

    std::string::size_type energy_pos = line.find("FINAL SINGLE");
    if (energy_pos != std::string::npos) {

      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      std::string energy = results[4];
      boost::trim(energy);
      orbitals.setQMEnergy(boost::lexical_cast<double>(energy));
      XTP_LOG(logDEBUG, *_pLog) << (boost::format("QM energy[Hrt]: %4.6f ") %
                                    orbitals.getDFTTotalEnergy())
                                       .str()
                                << flush;
    }

    std::string::size_type HFX_pos = line.find("Fraction HF Exchange ScalHFX");
    if (HFX_pos != std::string::npos) {
      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      double ScaHFX = boost::lexical_cast<double>(results.back());
      orbitals.setScaHFX(ScaHFX);
      XTP_LOG(logDEBUG, *_pLog)
          << "DFT with " << ScaHFX << " of HF exchange!" << flush;
    }

    std::string::size_type dim_pos = line.find("Basis Dimension");
    if (dim_pos != std::string::npos) {
      boost::algorithm::split(results, line, boost::is_any_of(" "),
                              boost::algorithm::token_compress_on);
      std::string dim =
          results[4];  // The 4th element of results vector is the Basis Dim
      boost::trim(dim);
      levels = boost::lexical_cast<int>(dim);
      XTP_LOG(logDEBUG, *_pLog) << "Basis Dimension: " << levels << flush;
      XTP_LOG(logDEBUG, *_pLog) << "Energy levels: " << levels << flush;
    }

    std::string::size_type OE_pos = line.find("ORBITAL ENERGIES");
    if (OE_pos != std::string::npos) {

      number_of_electrons = 0;
      getline(input_file, line);
      getline(input_file, line);
      getline(input_file, line);
      if (line.find("E(Eh)") == std::string::npos) {
        XTP_LOG(logDEBUG, *_pLog)
            << "Warning: Orbital Energies not found in log file" << flush;
      }
      for (int i = 0; i < levels; i++) {
        results = GetLineAndSplit(input_file, " ");
        std::string no = results[0];
        boost::trim(no);
        int levelnumber = boost::lexical_cast<int>(no);
        if (levelnumber != i) {
          XTP_LOG(logDEBUG, *_pLog) << "Have a look at the orbital energies "
                                       "something weird is going on"
                                    << flush;
        }
        std::string oc = results[1];
        boost::trim(oc);
        double occ = boost::lexical_cast<double>(oc);
        // We only count alpha electrons, each orbital must be empty or doubly
        // occupied
        if (occ == 2 || occ == 1) {
          number_of_electrons++;
          occupancy[i] = occ;
        } else if (occ == 0) {
          occupancy[i] = occ;
        } else {
          if (occ == 1) {
            XTP_LOG(logDEBUG, *_pLog)
                << "Watch out! No distinction between alpha and beta "
                   "electrons. Check if occ = 1 is suitable for your "
                   "calculation "
                << flush;
            number_of_electrons++;
            occupancy[i] = occ;
          } else {
            throw runtime_error(
                "Only empty or doubly occupied orbitals are allowed not "
                "running the right kind of DFT calculation");
          }
        }
        std::string e = results[2];
        boost::trim(e);
        energies[i] = boost::lexical_cast<double>(e);
      }
    }

    std::string::size_type success =
        line.find("*                     SUCCESS                       *");
    if (success != std::string::npos) {
      found_success = true;
    }
  }

  XTP_LOG(logDEBUG, *_pLog)
      << "Alpha electrons: " << number_of_electrons << flush;
  int occupied_levels = number_of_electrons;
  int unoccupied_levels = levels - occupied_levels;
  XTP_LOG(logDEBUG, *_pLog) << "Occupied levels: " << occupied_levels << flush;
  XTP_LOG(logDEBUG, *_pLog)
      << "Unoccupied levels: " << unoccupied_levels << flush;

  /************************************************************/

  // copying information to the orbitals object

  orbitals.setBasisSetSize(levels);
  orbitals.setNumberOfAlphaElectrons(number_of_electrons);
  orbitals.setNumberOfOccupiedLevels(occupied_levels);

  // copying energies to a vector
  orbitals.MOs().eigenvalues().resize(levels);
  //_level = 1;
  for (int i = 0; i < levels; i++) {
    orbitals.MOs().eigenvalues()[i] = energies[i];
  }

  XTP_LOG(logDEBUG, *_pLog) << "Done reading Log file" << flush;

  return found_success;
}
template <class T>
void Orca::GetCoordinates(T& mol, string& line, ifstream& input_file) const {
  std::string::size_type coordinates_pos =
      line.find("CARTESIAN COORDINATES (ANGSTROEM)");

  typedef typename std::iterator_traits<typename T::iterator>::value_type Atom;

  if (coordinates_pos != std::string::npos) {
    XTP_LOG(logDEBUG, *_pLog) << "Getting the coordinates" << flush;
    bool has_QMAtoms = mol.size() > 0;
    // three garbage lines
    getline(input_file, line);
    // now starts the data in format
    // _id type Qnuc x y z
    vector<string> row = GetLineAndSplit(input_file, "\t ");
    int nfields = row.size();
    int atom_id = 0;
    while (nfields == 4) {
      string atom_type = row.at(0);
      double x = boost::lexical_cast<double>(row.at(1));
      double y = boost::lexical_cast<double>(row.at(2));
      double z = boost::lexical_cast<double>(row.at(3));
      row = GetLineAndSplit(input_file, "\t ");
      nfields = row.size();
      Eigen::Vector3d pos(x, y, z);
      pos *= tools::conv::ang2bohr;
      if (has_QMAtoms == false) {
        mol.push_back(Atom(atom_id, atom_type, pos));
      } else {
        Atom& pAtom = mol.at(atom_id);
        pAtom.setPos(pos);
      }
      atom_id++;
    }
  }
}

bool Orca::CheckLogFile() {
  // check if the log file exists
  ifstream input_file(_run_dir + "/" + _log_file_name);

  if (input_file.fail()) {
    XTP_LOG(logERROR, *_pLog) << "Orca LOG is not found" << flush;
    return false;
  };

  std::string line;
  while (input_file) {
    getline(input_file, line);
    boost::trim(line);
    std::string::size_type error = line.find("FATAL ERROR ENCOUNTERED");
    if (error != std::string::npos) {
      XTP_LOG(logERROR, *_pLog) << "ORCA encountered a fatal error, maybe a "
                                   "look in the log file may help."
                                << flush;
      return false;
    }
    error = line.find(
        "mpirun detected that one or more processes exited with non-zero "
        "status");
    if (error != std::string::npos) {
      XTP_LOG(logERROR, *_pLog)
          << "ORCA had an mpi problem, maybe your openmpi version is not good."
          << flush;
      return false;
    }
  }
  return true;
}

// Parses the Orca gbw file and stores data in the Orbitals object

bool Orca::ParseMOsFile(Orbitals& orbitals) {
  if (!CheckLogFile()) return false;
  std::vector<double> coefficients;
  int basis_size = orbitals.getBasisSetSize();
  if (basis_size == 0) {
    throw runtime_error(
        "Basis size not set, calculator does not parse log file first");
  }

  XTP_LOG(logDEBUG, *_pLog)
      << "Reading the gbw file, this may or may not work so be careful: "
      << flush;
  ifstream infile;
  infile.open(_run_dir + "/" + _mo_file_name, ios::binary | ios::in);
  if (!infile) {
    throw runtime_error("Could not open " + _mo_file_name + " file");
  }
  infile.seekg(24, ios::beg);
  std::array<char, 8> buffer;
  infile.read(buffer.data(), 8);
  if (!infile) {
    infile.close();
    return false;
  }
  long int offset = *((long int*)buffer.data());

  infile.seekg(offset, ios::beg);
  infile.read(buffer.data(), 4);
  if (!infile) {
    infile.close();
    return false;
  }
  int op_read = *((int*)buffer.data());
  infile.seekg(offset + 4, ios::beg);
  infile.read(buffer.data(), 4);
  if (!infile) {
    infile.close();
    return false;
  }
  int dim_read = *((int*)buffer.data());
  infile.seekg(offset + 8, ios::beg);
  XTP_LOG(logDEBUG, *_pLog) << "Number of operators: " << op_read
                            << " Basis dimension: " << dim_read << flush;
  int n = op_read * dim_read * dim_read;
  for (int i = 0; i < n; i++) {
    infile.read(buffer.data(), 8);
    if (!infile) {
      infile.close();
      return false;
    }
    double mocoeff = *((double*)buffer.data());
    coefficients.push_back(mocoeff);
  }

  infile.close();
  // i -> MO, j -> AO
  orbitals.MOs().eigenvectors().resize(basis_size, basis_size);
  for (int i = 0; i < basis_size; i++) {
    for (int j = 0; j < basis_size; j++) {
      orbitals.MOs().eigenvectors()(j, i) = coefficients[j * basis_size + i];
    }
  }
  ReorderOutput(orbitals);
  XTP_LOG(logDEBUG, *_pLog) << "Done parsing" << flush;
  return true;
}

std::string Orca::getLName(int lnum) {
  if (lnum == 0) {
    return "S";
  } else if (lnum == 1) {
    return "P";
  } else if (lnum == 2) {
    return "D";
  } else if (lnum == 3) {
    return "F";
  } else {
    throw runtime_error(
        "Orca::getLName functions higher than F not implemented");
  }
  return "0";
}

std::string Orca::indent(const double& number) {
  std::stringstream ssnumber;
  if (number >= 0) {
    ssnumber << "    ";
  } else {
    ssnumber << "   ";
  }
  ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific
           << number;
  std::string snumber = ssnumber.str();
  return snumber;
}

}  // namespace xtp
}  // namespace votca
