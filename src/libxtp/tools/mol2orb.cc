#include "mol2orb.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmatom.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Mol2Orb::Initialize(const tools::Property& user_options) {
  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _job_name = options.ifExistsReturnElseReturnDefault<std::string>("job_name",
                                                                   _job_name);
  _moldenfile = _job_name + ".molden.input";
  _orbfile = _job_name + ".orb";
  _xyzfile = _job_name + ".xyz";

  _basisset_name = options.get(".basisset").as<std::string>();
  _aux_basisset_name = options.get(".auxbasisset").as<std::string>();
}

void Mol2Orb::addBasissetInfo(Orbitals& orbitals) {
  // We only need the basis set for checking correctness
  BasisSet bs;
  bs.Load(_basisset_name);

  // fill DFT AO basis by going through all atoms
  _basis.Fill(bs, orbitals.QMAtoms());
  orbitals.setDFTbasisName(_basisset_name);
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Basisset size: " << _basis.AOBasisSize()
      << std::flush;
  orbitals.setBasisSetSize(_basis.AOBasisSize());
  orbitals.setAuxbasisName(_aux_basisset_name);
}

// The returned string contains the line with the next section header
// or if it is the last section an empty string.
std::string Mol2Orb::readAtoms(QMMolecule& mol, const std::string& units,
                               std::ifstream& input_file) const {
  std::string line;
  std::istringstream iss(" ");
  while (std::getline(input_file, line)) {
    boost::trim(line);
    if (line == "" || line[0] == '[') {
      return line;
    }
    iss.str(line);
    iss.clear();

    // extract data
    double x, y, z;
    int atom_id;
    std::string junk;
    std::string atom_type;
    iss >> atom_type >> atom_id >> junk >> x >> y >> z;
    atom_id =
        atom_id - 1;  // molden uses indexing from 1, we use indexing from 0

    // Add data to orbitals object
    Eigen::Vector3d pos(x, y, z);
    if (units == "Angs") {
      pos = tools::conv::ang2bohr * pos;
    }
    mol.push_back(QMAtom(atom_id, atom_type, pos));
  }
  return "";
}

// The returned string contains the line with the next section header
// or, if it is the last section, an empty string.
std::string Mol2Orb::readMOs(Orbitals& orbitals, std::ifstream& input_file) {

  // setup space to store everything
  Index basis_size = orbitals.getBasisSetSize();
  int number_of_electrons = 0;
  if (basis_size == 0) {
    throw std::runtime_error(
        "Basis size not set, atoms were not parsed first.");
  }
  orbitals.MOs().eigenvalues().resize(basis_size);
  orbitals.MOs().eigenvectors().resize(basis_size, basis_size);

  // actual parsing
  std::string line;
  std::string tempStr;
  double tempDouble;
  Index tempIndex;
  std::istringstream iss(" ");
  for (int i = 0; i < basis_size; i++) {  // loop over mo's
    std::getline(input_file, line);       // skip symmetry label
    // energy line
    std::getline(input_file, line);
    iss.str(line);
    iss.clear();
    iss >> tempStr >> tempDouble;
    orbitals.MOs().eigenvalues()[i] = tempDouble;
    // spin channel line
    std::getline(input_file, line);
    iss.str(line);
    iss.clear();
    iss >> tempStr >> tempStr;
    if (tempStr == "Beta") {
      throw std::runtime_error(
          "Open shell systems are currently not supported");
    }
    // occupation line
    std::getline(input_file, line);
    iss.str(line);
    iss.clear();
    iss >> tempStr >> tempDouble;
    number_of_electrons += (int)tempDouble;

    // MO coefficients
    for (int j = 0; j < basis_size; j++) {  // loop over ao's
      std::getline(input_file, line);
      iss.str(line);
      iss.clear();
      iss >> tempIndex >> tempDouble;
      orbitals.MOs().eigenvectors()(j, i) = tempDouble;
    }
  }

  orbitals.setNumberOfAlphaElectrons(number_of_electrons);
  orbitals.setNumberOfOccupiedLevels(number_of_electrons / 2);

  // Get the correct ordering of atomic orbitals
  // Note that the molden specification specifies its own orbital ordering
  // Hence we can't use the qmpackage's ordering.
  OrbReorder reorder(_transpositions, _multipliers);
  reorder.reorderOrbitals(orbitals.MOs().eigenvectors(), _basis);

  getline(input_file, line);
  return line;
}

bool Mol2Orb::Evaluate() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;

  // Open the file
  std::ifstream input_file(_moldenfile);
  // Check if succesfull
  if (input_file.fail()) {
    XTP_LOG(Log::error, _log) << "Could not open " << _moldenfile << std::flush;
    return false;
  }

  XTP_LOG(Log::error, _log)
      << "Started parsing " << _moldenfile << " to " << _orbfile << std::flush;

  QMMolecule& mol = orbitals.QMAtoms();

  std::string line;
  std::getline(input_file, line);
  while (input_file) {
    boost::trim(line);
    if (line[0] != '[') {
      // ignore non-relevant lines
      std::getline(input_file, line);
      continue;
    }

    // Extract the part between square brackets
    long unsigned close = line.find("]");
    std::string sectionType = line.substr(1, close - 1);

    // Import data from relevant sections
    if (sectionType == "Atoms") {
      std::string units = line.substr(close + 1);  // extra units Atomic Units
                                                   // (AU) or Angstrom (Angs)
      boost::trim(units);
      XTP_LOG(Log::error, _log)
          << "Reading atoms using " << units << " units." << std::flush;
      line = readAtoms(mol, units, input_file);
      addBasissetInfo(orbitals);
    } else if (sectionType == "GTO") {
      XTP_LOG(Log::error, _log)
          << "Basisset specification is ignored." << std::flush;
      XTP_LOG(Log::error, _log)
          << "Basissets are specified via the mol2orb.xml options file."
          << std::flush;
      std::getline(input_file, line);
    } else if (sectionType == "MO") {
      if (mol.size() == 0) {
        XTP_LOG(Log::error, _log)
            << "Atoms should be specified before MO coefficients!"
            << std::flush;
        return false;
      } else {
        XTP_LOG(Log::error, _log)
            << "Reading molecular orbital coefficients" << std::flush;
        line = readMOs(orbitals, input_file);
      }
    } else if (sectionType == "STO") {
      XTP_LOG(Log::error, _log)
          << "Slater Type Orbitals (STOs) are not supported in VOTCA-XTP."
          << std::flush;
      std::getline(input_file, line);
    } else {
      XTP_LOG(Log::error, _log)
          << "VOTCA-XTP ignores " << sectionType << " section." << std::flush;
      std::getline(input_file, line);
    }
  }

  // Save orbitals object
  XTP_LOG(Log::error, _log) << "Saving data to " << _orbfile << std::flush;
  orbitals.WriteToCpt(_orbfile);
  XTP_LOG(Log::error, _log) << "Done parsing! \n" << std::flush;
  return true;
}

}  // namespace xtp
}  // namespace votca