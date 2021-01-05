#include <sstream>

#include "votca/xtp/basisset.h"
#include "votca/xtp/molden.h"
#include <boost/algorithm/string.hpp>

namespace votca {
namespace xtp {

void Molden::writeAtoms(const Orbitals& orbitals,
                        std::ofstream& outFile) const {
  for (const auto& atom : orbitals.QMAtoms()) {
    const Eigen::Vector3d& pos = atom.getPos();
    outFile << boost::format("%-4s %5d %5d %22.12e %22.10e %22.10e\n") %
                   atom.getElement() % (atom.getId() + 1) %
                   atom.getElementNumber() % pos[0] % pos[1] % pos[2];
  }
}

void Molden::writeMOs(const Orbitals& orbitals, std::ofstream& outFile) const {
  Eigen::VectorXd energies = orbitals.MOs().eigenvalues();
  bool fromVotcaToExternal = true;
  OrbReorder reorder(_reorderList, _multipliers, fromVotcaToExternal);
  Eigen::MatrixXd moCoefficients = orbitals.MOs().eigenvectors();
  reorder.reorderOrbitals(moCoefficients, orbitals.SetupDftBasis());

  for (Index i = 0; i < orbitals.getBasisSetSize(); i++) {  // over columns
    outFile << "Sym= \n";
    outFile << boost::format("Ene= %-20.12e\n") % energies[i];
    outFile << "Spin= Alpha\n";
    outFile << boost::format("Occup= %-5.2f\n") %
                   (2 * (i < orbitals.getLumo()));
    for (Index j = 0; j < orbitals.getBasisSetSize(); j++) {
      outFile << boost::format("%5d %22.12e\n") % (j + 1) %
                     moCoefficients(j, i);
    }
  }
}

void Molden::writeBasisSet(const Orbitals& orbitals,
                           std::ofstream& outFile) const {
  AOBasis basis = orbitals.SetupDftBasis();
  for (const auto& atom : orbitals.QMAtoms()) {
    // The 0 in the format string of the next line is meaningless it
    // is included for backwards compatibility of molden files
    outFile << boost::format("%4d 0 \n") % (atom.getId() + 1);
    const std::vector<const AOShell*> shells =
        basis.getShellsofAtom(atom.getId());
    for (const AOShell* shell : shells) {
      // The 1.0 at the end of the next line is meaningless it
      // is included for backwards compatibility of molden files
      outFile << boost::format("%-3s %4d %3.1f \n") %
                     boost::to_lower_copy(EnumToString(shell->getL())) %
                     shell->getSize() % 1.0;
      for (const AOGaussianPrimitive& gaussian : *shell) {
        outFile << boost::format("%22.10e %22.10e\n") % gaussian.getDecay() %
                       gaussian.getContraction();
      }
    }
    outFile << " \n";
  }
}

void Molden::WriteFile(const std::string& filename,
                       const Orbitals& orbitals) const {
  if (!orbitals.hasDFTbasisName()) {
    throw std::runtime_error(".orb file does not contain a basisset name");
  }

  std::ofstream outFile(filename);
  if (outFile.is_open()) {

    XTP_LOG(Log::error, _log) << "Writing data to " << filename << std::flush;

    // print Header
    outFile << "[Molden Format]\n";
    outFile << "[Title]\n";
    outFile << "Molden file created by VOTCA-XTP\n";
    outFile << " \n";

    outFile << "[Atoms] AU\n";
    writeAtoms(orbitals, outFile);

    outFile << "[GTO] \n";
    writeBasisSet(orbitals, outFile);

    // indicate spherical D F and G functions
    outFile << "[5D] \n[7F] \n[9G] \n";

    outFile << "[MO]\n";
    writeMOs(orbitals, outFile);
    XTP_LOG(Log::error, _log)
        << "Finished writing to molden file." << std::flush;
  }
}

std::string Molden::readAtoms(QMMolecule& mol, const std::string& units,
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
    Index atom_id;
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
std::string Molden::readMOs(Orbitals& orbitals,
                            std::ifstream& input_file) const {

  // setup space to store everything
  Index basis_size = orbitals.getBasisSetSize();
  Index number_of_electrons = 0;
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
  for (Index i = 0; i < basis_size; i++) {  // loop over mo's
    std::getline(input_file, line);         // skip symmetry label
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

  OrbReorder reorder(_reorderList, _multipliers);
  reorder.reorderOrbitals(orbitals.MOs().eigenvectors(),
                          orbitals.SetupDftBasis());

  getline(input_file, line);
  return line;
}

void Molden::addBasissetInfo(Orbitals& orbitals) const {
  orbitals.setDFTbasisName(_basisset_name);
  orbitals.setBasisSetSize(orbitals.SetupDftBasis().AOBasisSize());
  orbitals.setAuxbasisName(_aux_basisset_name);
}

void Molden::parseMoldenFile(const std::string& filename,
                             Orbitals& orbitals) const {

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