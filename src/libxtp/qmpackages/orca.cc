
/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Standard includes
#include <cstdio>
#include <iomanip>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <stdexcept>
#include <string>
#include <votca/tools/elements.h>
#include <votca/tools/getline.h>

// Local VOTCA includes
#include "votca/tools/globals.h"
#include "votca/tools/tokenizer.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/molden.h"
#include "votca/xtp/orbitals.h"

// Local private VOTCA includes
#include "orca.h"

namespace votca {
namespace xtp {
using namespace std;

void Orca::ParseSpecificOptions(const tools::Property& options) {

  // Orca file names
  std::string fileName = options.get("temporary_file").as<std::string>();

  input_file_name_ = fileName + ".inp";
  log_file_name_ = fileName + ".log";
  shell_file_name_ = fileName + ".sh";
  mo_file_name_ = fileName + ".gbw";
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
  XTP_LOG(Log::error, *pLog_) << "Loaded Basis Set " << bs_name << flush;
  ofstream el_file;

  el_file.open(el_file_name);
  el_file << "$DATA" << endl;

  for (const std::string& element_name : UniqueElements) {
    const Element& element = bs.getElement(element_name);
    el_file << elementInfo.getEleFull(element_name) << endl;
    for (const Shell& shell : element) {
      el_file << xtp::EnumToString(shell.getL()) << " " << shell.getSize()
              << endl;
      Index sh_idx = 0;
      for (const GaussianPrimitive& gaussian : shell) {
        sh_idx++;
        el_file << " " << sh_idx << " " << indent(gaussian.decay());
        el_file << " " << indent(gaussian.contraction());
        el_file << endl;
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
             << setiosflags(ios::fixed) << setprecision(6) << pos.x()
             << setw(12) << setiosflags(ios::fixed) << setprecision(6)
             << pos.y() << setw(12) << setiosflags(ios::fixed)
             << setprecision(6) << pos.z() << endl;
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
  ecp.Load(options_.get("ecp").as<std::string>());

  XTP_LOG(Log::error, *pLog_) << "Loaded Pseudopotentials "
                              << options_.get("ecp").as<std::string>() << flush;

  for (const std::string& element_name : UniqueElements) {
    try {
      ecp.getElement(element_name);
    } catch (std::runtime_error& error) {
      XTP_LOG(Log::error, *pLog_)
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
             << " " << EnumToString(element.getLmax()) << endl;
    // For Orca the order doesn't matter but let's write it in ascending order
    // write remaining shells in ascending order s,p,d...
    for (Index i = 0; i <= Index(element.getLmax()); i++) {
      for (const ECPShell& shell : element) {
        if (Index(shell.getL()) == i) {
          // shell type, number primitives, scale factor
          inp_file << xtp::EnumToString(shell.getL()) << " " << shell.getSize()
                   << endl;
          Index sh_idx = 0;
          for (const ECPGaussianPrimitive& gaussian : shell) {
            sh_idx++;
            inp_file << sh_idx << " " << gaussian.decay_ << " "
                     << gaussian.contraction_ << " " << gaussian.power_ << endl;
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
  this->options_.addTree("orca.pointcharges", "\"background.crg\"");
}

/* For QM/MM the molecules in the MM environment are represented by
 * their atomic partial charge distributions. ORCA expects them in
 * q,x,y,z format in a separate file "background.crg"
 */
void Orca::WriteBackgroundCharges() {

  std::ofstream crg_file;
  std::string crg_file_name_full_ = run_dir_ + "/background.crg";
  crg_file.open(crg_file_name_full_);
  Index total_background = 0;

  for (const std::unique_ptr<StaticSite>& site : externalsites_) {
    if (site->getCharge() != 0.0) {
      total_background++;
    }
    total_background += SplitMultipoles(*site).size();
  }  // counting only

  crg_file << total_background << endl;
  boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
  // now write
  for (const std::unique_ptr<StaticSite>& site : externalsites_) {
    Eigen::Vector3d pos = site->getPos() * tools::conv::bohr2ang;
    string sitestring =
        boost::str(fmt % site->getCharge() % pos.x() % pos.y() % pos.z());
    if (site->getCharge() != 0.0) {
      crg_file << sitestring << endl;
    }
    std::vector<MinimalMMCharge> split_multipoles = SplitMultipoles(*site);
    for (const auto& mpoles : split_multipoles) {
      Eigen::Vector3d pos2 = mpoles.pos_ * tools::conv::bohr2ang;
      string multipole =
          boost::str(fmt % mpoles.q_ % pos2.x() % pos2.y() % pos2.z());
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
  std::string scratch_dir_backup = scratch_dir_;
  std::ofstream inp_file;
  std::string inp_file_name_full = run_dir_ + "/" + input_file_name_;
  inp_file.open(inp_file_name_full);
  // header
  inp_file << "* xyz  " << charge_ << " " << spin_ << endl;
  Index threads = OPENMP::getMaxThreads();
  const QMMolecule& qmatoms = orbitals.QMAtoms();
  // put coordinates
  WriteCoordinates(inp_file, qmatoms);
  // add parallelization info
  inp_file << "%pal\n"
           << "nprocs " << threads << "\nend"
           << "\n"
           << endl;
  // basis set info
  std::string el_file_name = run_dir_ + "/" + "system.bas";
  WriteBasisset(qmatoms, basisset_name_, el_file_name);
  inp_file << "%basis\n";
  inp_file << "GTOName"
           << " "
           << "="
           << "\"system.bas\";" << endl;
  if (options_.exists("auxbasisset")) {
    std::string aux_file_name = run_dir_ + "/" + "system.aux";
    std::string auxbasisset_name =
        options_.get("auxbasisset").as<std::string>();
    WriteBasisset(qmatoms, auxbasisset_name, aux_file_name);
    inp_file << "GTOAuxName"
             << " "
             << "="
             << "\"system.aux\";" << endl;
  }  // write_auxbasis set

  // ECPs
  if (options_.exists("ecp")) {
    WriteECP(inp_file, qmatoms);
  }
  inp_file << "end\n "
           << "\n"
           << endl;  // This end is for the basis set block
  if (!externalsites_.empty()) {
    WriteBackgroundCharges();
  }

  // External Electric field
  if (options_.exists("externalfield")) {
    Eigen::Vector3d field = options_.get("externalfield").as<Eigen::Vector3d>();
    inp_file << "%scf\n ";
    inp_file << "  efield " << field.x() << ", " << field.y() << ", "
             << field.z() << "\n";
    inp_file << "end\n";
    inp_file << std::endl;
  }
  std::string input_options;
  // Write Orca section specified by the user
  for (const auto& prop : this->options_.get("orca")) {
    const std::string& prop_name = prop.name();
    if (prop_name == "pointcharges") {
      input_options += this->CreateInputSection("orca.pointcharges");
    } else if (prop_name != "method") {
      input_options += this->CreateInputSection("orca." + prop_name);
    }
  }
  // Write main DFT method
  input_options += this->WriteMethod();
  inp_file << input_options;
  inp_file << '\n';
  inp_file.close();
  // and now generate a shell script to run both jobs, if neccessary

  XTP_LOG(Log::info, *pLog_)
      << "Setting the scratch dir to " << scratch_dir_ + temp_suffix << flush;
  scratch_dir_ = scratch_dir_backup + temp_suffix;
  WriteShellScript();
  scratch_dir_ = scratch_dir_backup;
  return true;
}

bool Orca::WriteShellScript() {
  ofstream shell_file;
  std::string shell_file_name_full = run_dir_ + "/" + shell_file_name_;
  shell_file.open(shell_file_name_full);
  shell_file << "#!/bin/bash" << endl;
  shell_file << "mkdir -p " << scratch_dir_ << endl;

  if (options_.get("initial_guess").as<std::string>() == "orbfile") {
    if (!(boost::filesystem::exists(run_dir_ + "/molA.gbw") &&
          boost::filesystem::exists(run_dir_ + "/molB.gbw"))) {
      throw runtime_error(
          "Using guess relies on a molA.gbw and a molB.gbw file being in the "
          "directory.");
    }
    shell_file << options_.get("executable").as<std::string>()
               << "_mergefrag molA.gbw molB.gbw dimer.gbw > merge.log" << endl;
  }
  shell_file << options_.get("executable").as<std::string>() << " "
             << input_file_name_ << " > " << log_file_name_
             << endl;  //" 2> run.error" << endl;
  std::string base_name = mo_file_name_.substr(0, mo_file_name_.size() - 4);
  shell_file << options_.get("executable").as<std::string>() << "_2mkl "
             << base_name << " -molden" << endl;
  shell_file.close();
  return true;
}

/**
 * Runs the Orca job.
 */
bool Orca::RunDFT() {

  XTP_LOG(Log::error, *pLog_) << "Running Orca job\n" << flush;

  if (std::system(nullptr)) {

    std::string command = "cd " + run_dir_ + "; sh " + shell_file_name_;
    Index check = std::system(command.c_str());
    if (check == -1) {
      XTP_LOG(Log::error, *pLog_)
          << input_file_name_ << " failed to start" << flush;
      return false;
    }
    if (CheckLogFile()) {
      XTP_LOG(Log::error, *pLog_) << "Finished Orca job" << flush;
      return true;
    } else {
      XTP_LOG(Log::error, *pLog_) << "Orca job failed" << flush;
    }
  } else {
    XTP_LOG(Log::error, *pLog_)
        << input_file_name_ << " failed to start" << flush;
    return false;
  }

  return true;
}

/**
 * Cleans up after the Orca job
 */
void Orca::CleanUp() {

  if (options_.get("initial_guess").as<std::string>() == "orbfile") {
    remove((run_dir_ + "/" + "molA.gbw").c_str());
    remove((run_dir_ + "/" + "molB.gbw").c_str());
    remove((run_dir_ + "/" + "dimer.gbw").c_str());
  }
  // cleaning up the generated files
  if (!cleanup_.empty()) {
    tools::Tokenizer tok_cleanup(cleanup_, ",");
    for (const std::string& substring : tok_cleanup) {
      if (substring == "inp") {
        std::string file_name = run_dir_ + "/" + input_file_name_;
        remove(file_name.c_str());
      }

      if (substring == "bas") {
        std::string file_name = run_dir_ + "/system.bas";
        remove(file_name.c_str());
      }

      if (substring == "log") {
        std::string file_name = run_dir_ + "/" + log_file_name_;
        remove(file_name.c_str());
      }

      if (substring == "gbw") {
        std::string file_name = run_dir_ + "/" + mo_file_name_;
        remove(file_name.c_str());
      }

      if (substring == "ges") {
        std::string file_name = run_dir_ + "/system.ges";
        remove(file_name.c_str());
      }
      if (substring == "prop") {
        std::string file_name = run_dir_ + "/system.prop";
        remove(file_name.c_str());
      }
    }
  }
  return;
}

StaticSegment Orca::GetCharges() const {

  StaticSegment result("charges", 0);

  XTP_LOG(Log::error, *pLog_) << "Parsing " << log_file_name_ << flush;
  std::string log_file_name_full = run_dir_ + "/" + log_file_name_;
  std::string line;

  std::ifstream input_file(log_file_name_full);
  while (input_file) {
    tools::getline(input_file, line);
    boost::trim(line);
    GetCoordinates(result, line, input_file);

    std::string::size_type charge_pos = line.find("CHELPG Charges");

    if (charge_pos != std::string::npos) {
      XTP_LOG(Log::error, *pLog_) << "Getting charges" << flush;
      tools::getline(input_file, line);
      std::vector<std::string> row = GetLineAndSplit(input_file, "\t ");
      Index nfields = Index(row.size());
      bool hasAtoms = result.size() > 0;
      while (nfields == 4) {
        Index atom_id = boost::lexical_cast<Index>(row.at(0));
        std::string atom_type = row.at(1);
        double atom_charge = boost::lexical_cast<double>(row.at(3));
        row = GetLineAndSplit(input_file, "\t ");
        nfields = Index(row.size());
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
  ifstream input_file((run_dir_ + "/" + log_file_name_));
  bool has_pol = false;

  Eigen::Matrix3d pol = Eigen::Matrix3d::Zero();
  while (input_file) {
    tools::getline(input_file, line);
    boost::trim(line);

    std::string::size_type pol_pos = line.find("THE POLARIZABILITY TENSOR");
    if (pol_pos != std::string::npos) {
      XTP_LOG(Log::error, *pLog_) << "Getting polarizability" << flush;
      tools::getline(input_file, line);
      tools::getline(input_file, line);
      tools::getline(input_file, line);

      if (line.find("The raw cartesian tensor (atomic units)") ==
          std::string::npos) {
        throw std::runtime_error(
            "Could not find cartesian polarization tensor");
      }

      for (Index i = 0; i < 3; i++) {
        tools::getline(input_file, line);
        std::vector<double> values =
            tools::Tokenizer(line, " ").ToVector<double>();
        if (values.size() != 3) {
          throw std::runtime_error("polarization line " + line +
                                   " cannot be parsed");
        }
        Eigen::Vector3d row;
        row << values[0], values[1], values[2];
        pol.row(i) = row;
      }

      has_pol = true;
    }
  }
  if (!has_pol) {
    throw std::runtime_error("Could not find polarization in logfile");
  }
  return pol;
}

bool Orca::ParseLogFile(Orbitals& orbitals) {
  bool found_success = false;
  orbitals.setQMpackage(getPackageName());
  orbitals.setDFTbasisName(basisset_name_);
  if (options_.exists("ecp")) {
    orbitals.setECPName(options_.get("ecp").as<std::string>());
  }

  orbitals.setXCGrid("xfine");  // TODO find a better approximation for the
                                // orca grid.
  orbitals.setXCFunctionalName(options_.get("functional").as<std::string>());

  XTP_LOG(Log::error, *pLog_) << "Parsing " << log_file_name_ << flush;
  std::string log_file_name_full = run_dir_ + "/" + log_file_name_;
  // check if LOG file is complete
  if (!CheckLogFile()) {
    return false;
  }
  std::map<Index, double> energies;
  std::map<Index, double> occupancy;

  std::string line;
  Index levels = 0;
  Index number_of_electrons = 0;
  std::vector<std::string> results;

  std::ifstream input_file(log_file_name_full);

  if (input_file.fail()) {
    XTP_LOG(Log::error, *pLog_)
        << "File " << log_file_name_full << " not found " << flush;
    return false;
  } else {
    XTP_LOG(Log::error, *pLog_)
        << "Reading Coordinates and occupationnumbers and energies from "
        << log_file_name_full << flush;
  }
  // Coordinates of the final configuration depending on whether it is an
  // optimization or not

  QMMolecule& mol = orbitals.QMAtoms();
  while (input_file) {
    tools::getline(input_file, line);
    boost::trim(line);

    GetCoordinates(mol, line, input_file);

    std::string::size_type energy_pos = line.find("FINAL SINGLE");
    if (energy_pos != std::string::npos) {
      results = tools::Tokenizer(line, " ").ToVector();
      std::string energy = results[4];
      boost::trim(energy);
      orbitals.setQMEnergy(boost::lexical_cast<double>(energy));
      XTP_LOG(Log::error, *pLog_) << (boost::format("QM energy[Hrt]: %4.6f ") %
                                      orbitals.getDFTTotalEnergy())
                                         .str()
                                  << flush;
    }

    std::string::size_type HFX_pos = line.find("Fraction HF Exchange ScalHFX");

    if (HFX_pos != std::string::npos) {
      results = tools::Tokenizer(line, " ").ToVector();
      double ScaHFX = boost::lexical_cast<double>(results.back());

      orbitals.setScaHFX(ScaHFX);
      XTP_LOG(Log::error, *pLog_)
          << "DFT with " << ScaHFX << " of HF exchange!" << flush;
    }

    std::string::size_type dim_pos = line.find("Basis Dimension");
    if (dim_pos != std::string::npos) {
      results = tools::Tokenizer(line, " ").ToVector();
      std::string dim =
          results[4];  // The 4th element of results vector is the Basis Dim
      boost::trim(dim);
      levels = boost::lexical_cast<Index>(dim);
      XTP_LOG(Log::info, *pLog_) << "Basis Dimension: " << levels << flush;
      XTP_LOG(Log::info, *pLog_) << "Energy levels: " << levels << flush;
    }

    std::string::size_type OE_pos = line.find("ORBITAL ENERGIES");
    if (OE_pos != std::string::npos) {

      number_of_electrons = 0;
      tools::getline(input_file, line);
      tools::getline(input_file, line);
      tools::getline(input_file, line);
      if (line.find("E(Eh)") == std::string::npos) {
        XTP_LOG(Log::error, *pLog_)
            << "Warning: Orbital Energies not found in log file" << flush;
      }
      for (Index i = 0; i < levels; i++) {
        results = GetLineAndSplit(input_file, " ");
        std::string no = results[0];
        boost::trim(no);
        Index levelnumber = boost::lexical_cast<Index>(no);
        if (levelnumber != i) {
          XTP_LOG(Log::error, *pLog_) << "Have a look at the orbital energies "
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
          if (occ == 1) {
            XTP_LOG(Log::error, *pLog_)
                << "Watch out! No distinction between alpha and beta "
                   "electrons. Check if occ = 1 is suitable for your "
                   "calculation "
                << flush;
          }
        } else if (occ == 0) {
          occupancy[i] = occ;
        } else {
          throw runtime_error(
              "Only empty or doubly occupied orbitals are allowed not "
              "running the right kind of DFT calculation");
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

  XTP_LOG(Log::info, *pLog_)
      << "Alpha electrons: " << number_of_electrons << flush;
  Index occupied_levels = number_of_electrons;
  Index unoccupied_levels = levels - occupied_levels;
  XTP_LOG(Log::info, *pLog_) << "Occupied levels: " << occupied_levels << flush;
  XTP_LOG(Log::info, *pLog_)
      << "Unoccupied levels: " << unoccupied_levels << flush;

  /************************************************************/

  // copying information to the orbitals object

  orbitals.setBasisSetSize(levels);
  orbitals.setNumberOfAlphaElectrons(number_of_electrons);
  orbitals.setNumberOfOccupiedLevels(occupied_levels);

  // copying energies to a vector
  orbitals.MOs().eigenvalues().resize(levels);
  // level_ = 1;
  for (Index i = 0; i < levels; i++) {
    orbitals.MOs().eigenvalues()[i] = energies[i];
  }

  XTP_LOG(Log::error, *pLog_) << "Done reading Log file" << flush;

  return found_success;
}
template <class T>
void Orca::GetCoordinates(T& mol, string& line, ifstream& input_file) const {
  std::string::size_type coordinates_pos =
      line.find("CARTESIAN COORDINATES (ANGSTROEM)");

  using Atom = typename T::Atom_Type;

  if (coordinates_pos != std::string::npos) {
    XTP_LOG(Log::error, *pLog_) << "Getting the coordinates" << flush;
    bool has_QMAtoms = mol.size() > 0;
    // three garbage lines
    tools::getline(input_file, line);
    // now starts the data in format
    //  id_ type Qnuc x y z
    vector<string> row = GetLineAndSplit(input_file, "\t ");
    Index nfields = Index(row.size());
    Index atom_id = 0;
    while (nfields == 4) {
      string atom_type = row.at(0);
      double x = boost::lexical_cast<double>(row.at(1));
      double y = boost::lexical_cast<double>(row.at(2));
      double z = boost::lexical_cast<double>(row.at(3));
      row = GetLineAndSplit(input_file, "\t ");
      nfields = Index(row.size());
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
  ifstream input_file(run_dir_ + "/" + log_file_name_);

  if (input_file.fail()) {
    XTP_LOG(Log::error, *pLog_) << "Orca LOG is not found" << flush;
    return false;
  };

  std::string line;
  while (input_file) {
    tools::getline(input_file, line);
    boost::trim(line);
    std::string::size_type error = line.find("FATAL ERROR ENCOUNTERED");
    if (error != std::string::npos) {
      XTP_LOG(Log::error, *pLog_) << "ORCA encountered a fatal error, maybe a "
                                     "look in the log file may help."
                                  << flush;
      return false;
    }
    error = line.find(
        "mpirun detected that one or more processes exited with non-zero "
        "status");
    if (error != std::string::npos) {
      XTP_LOG(Log::error, *pLog_)
          << "ORCA had an mpi problem, maybe your openmpi version is not good."
          << flush;
      return false;
    }
  }
  return true;
}

// Parses the molden file from orca and stores data in the Orbitals object
bool Orca::ParseMOsFile(Orbitals& orbitals) {
  if (!CheckLogFile()) {
    return false;
  }
  std::vector<double> coefficients;
  Index basis_size = orbitals.getBasisSetSize();
  if (basis_size == 0) {
    throw runtime_error(
        "Basis size not set, calculator does not parse log file first");
  }

  XTP_LOG(Log::error, *pLog_) << "Reading Molden file" << flush;

  Molden molden(*pLog_);

  if (orbitals.getDFTbasisName() == "") {
    throw runtime_error(
        "Basisset names should be set before reading the molden file.");
  }

  molden.setBasissetInfo(orbitals.getDFTbasisName());
  std::string file_name = run_dir_ + "/" +
                          mo_file_name_.substr(0, mo_file_name_.size() - 4) +
                          ".molden.input";
  XTP_LOG(Log::error, *pLog_) << "Molden file: " << file_name << flush;
  std::ifstream molden_file(file_name);
  if (!molden_file.good()) {
    throw std::runtime_error(
        "Could not find the molden input file for the MO coefficients.\nIf you "
        "have run the orca calculation manually or use data from an old\n"
        "calculation, make sure that besides the .gbw file a .molden.input is\n"
        "present. If not, convert the .gbw file to a .molden.input file with\n"
        "the orca_2mkl tool from orca.\nAn example, if you have a benzene.gbw "
        "file run:\n    orca_2mkl benzene -molden\n");
  }

  molden.parseMoldenFile(file_name, orbitals);

  XTP_LOG(Log::error, *pLog_) << "Done parsing" << flush;
  return true;
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

std::string Orca::CreateInputSection(const std::string& key) const {
  std::stringstream stream;
  std::string section = key.substr(key.find(".") + 1);
  stream << "%" << section;
  if (KeywordIsSingleLine(key)) {
    stream << " " << options_.get(key).as<std::string>() << "\n";
  } else {
    stream << "\n"
           << options_.get(key).as<std::string>() << "\n"
           << "end\n";
  }

  return stream.str();
}

bool Orca::KeywordIsSingleLine(const std::string& key) const {
  tools::Tokenizer values(this->options_.get(key).as<std::string>(), " ");
  std::vector<std::string> words = values.ToVector();
  return ((words.size() <= 1) ? true : false);
}

std::string Orca::WriteMethod() const {
  std::stringstream stream;
  std::string opt = (options_.get("optimize").as<bool>()) ? "Opt" : "";
  const tools::Property& orca = options_.get("orca");
  std::string user_method =
      (orca.exists("method")) ? orca.get("method").as<std::string>() : "";
  std::string convergence = "";
  if (!orca.exists("scf")) {
    convergence = this->convergence_map_.at(
        options_.get("convergence_tightness").as<std::string>());
  }
  stream << "! DFT " << this->GetOrcaFunctionalName() << " " << convergence
         << " " << opt
         << " "
         // additional properties provided by the user
         << user_method << "\n";
  return stream.str();
}

std::string Orca::GetOrcaFunctionalName() const {

  std::map<std::string, std::string> votca_to_orca;

  votca_to_orca["XC_HYB_GGA_XC_B1LYP"] = "B1LYP";
  votca_to_orca["XC_HYB_GGA_XC_B3LYP"] = "B3LYP";
  votca_to_orca["XC_HYB_GGA_XC_PBEH"] = "PBE0";
  votca_to_orca["XC_GGA_C_PBE"] = "PBE";
  votca_to_orca["XC_GGA_X_PBE"] = "PBE";

  std::string votca_functional =
      options_.get("functional").as<std::vector<std::string>>()[0];

  std::string orca_functional;
  if (votca_to_orca.count(votca_functional)) {
    orca_functional = votca_to_orca.at(votca_functional);
  } else if (options_.exists("orca." + votca_functional)) {
    orca_functional =
        options_.get("orca." + votca_functional).as<std::string>();
  } else {
    throw std::runtime_error(
        "Cannot translate " + votca_functional +
        " to orca functional names. Please add a <" + votca_functional +
        "> to your orca input options with the functional orca should use.");
  }
  return orca_functional;
}

}  // namespace xtp
}  // namespace votca
