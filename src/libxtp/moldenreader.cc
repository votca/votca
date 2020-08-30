#include <sstream>

#include "votca/xtp/moldenreader.h"

namespace votca {
namespace xtp {

// The returned string contains the line with the next section header
// or if it is the last section an empty string.
std::string MoldenReader::readAtoms(QMMolecule& mol, const std::string& units,
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
std::string MoldenReader::readMOs(Orbitals& orbitals,
                                  std::ifstream& input_file) const {

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

  OrbReorder reorder(_reorderList, _multipliers, true);
  reorder.reorderOrbitals(orbitals.MOs().eigenvectors(), _basis);

  getline(input_file, line);
  return line;
}

void MoldenReader::addBasissetInfo(Orbitals& orbitals) {
  BasisSet bs;
  bs.Load(_basisset_name);

  _basis.Fill(bs, orbitals.QMAtoms());
  orbitals.setDFTbasisName(_basisset_name);
  orbitals.setBasisSetSize(_basis.AOBasisSize());
  orbitals.setAuxbasisName(_aux_basisset_name);
}

void MoldenReader::parseMoldenFile(const std::string& filename,
                                   Orbitals& orbitals) {

  if (_basisset_name == "") {
    throw std::runtime_error(
        "Basisset names should be set before reading the molden file.");
  }

  std::ifstream input_file(filename);
  // Check if succesfull
  if (input_file.fail()) {
    throw std::runtime_error("Could not open molden file.");
  }

  std::string line;
  std::getline(input_file, line);
  while (input_file) {
    boost::trim(line);
    if (line[0] != '[') {  // ignore non-relevant lines
      std::getline(input_file, line);
      continue;
    }

    // Extract the part between square brackets
    long unsigned close = line.find("]");
    std::string sectionType = line.substr(1, close - 1);

    // Import data from relevant sections
    if (sectionType == "Atoms") {
      std::string units = line.substr(close + 1);
      boost::trim(units);
      XTP_LOG(Log::error, _log)
          << "Reading atoms using " << units << " units." << std::flush;
      line = readAtoms(orbitals.QMAtoms(), units, input_file);
      addBasissetInfo(orbitals);
    } else if (sectionType == "GTO") {
      XTP_LOG(Log::error, _log)
          << "Basisset specification is ignored." << std::flush;
      XTP_LOG(Log::error, _log)
          << "Basissets are specified via the mol2orb.xml options file."
          << std::flush;
      std::getline(input_file, line);
    } else if (sectionType == "MO") {
      if (orbitals.QMAtoms().size() == 0) {
        throw std::runtime_error(
            "Atoms should be specified before MO coefficients.");
      } else {
        XTP_LOG(Log::error, _log)
            << "Reading molecular orbital coefficients" << std::flush;
        line = readMOs(orbitals, input_file);
      }
    } else if (sectionType == "STO") {
      throw std::runtime_error(
          "Slater Type Orbitals (STOs) are not supported in VOTCA-XTP.");
    } else {
      std::getline(input_file, line);
    }
  }

  XTP_LOG(Log::error, _log) << "Done parsing molden file" << std::flush;
}

}  // namespace xtp
}  // namespace votca
