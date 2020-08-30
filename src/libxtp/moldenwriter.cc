#include "votca/xtp/moldenwriter.h"

#include "votca/xtp/basisset.h"
#include <boost/algorithm/string.hpp>

namespace votca {
namespace xtp {

void MoldenWriter::writeAtoms(const Orbitals& orbitals,
                              std::ofstream& outFile) const {
  for (const auto& atom : orbitals.QMAtoms()) {
    const Eigen::Vector3d& pos = atom.getPos();
    outFile << boost::format("%-4s %5d %5d %22.12e %22.10e %22.10e\n") %
                   atom.getElement() % (atom.getId() + 1) %
                   atom.getTotalNuccharge() % pos[0] % pos[1] % pos[2];
  }
}

void MoldenWriter::writeMOs(const Orbitals& orbitals,
                            std::ofstream& outFile) const {
  Eigen::VectorXd energies = orbitals.MOs().eigenvalues();
  OrbReorder reorder(_reorderList, _multipliers);
  Eigen::MatrixXd moCoefficients = orbitals.MOs().eigenvectors();
  reorder.reorderOrbitals(moCoefficients, _basis);

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

void MoldenWriter::writeBasisSet(const Orbitals& orbitals,
                                 std::ofstream& outFile) const {

  for (const auto& atom : orbitals.QMAtoms()) {
    const Element& element = _bs.getElement(atom.getElement());
    // The 0 in the format string of the next line is meaningless it
    // is included for backwards compatibility of molden files
    outFile << boost::format("%4d 0 \n") % (atom.getId() + 1);
    for (const Shell& shell : element) {
      // The 1.0 at the end of the next line is meaningless it
      // is included for backwards compatibility of molden files
      outFile << boost::format("%-3s %4d %3.1f \n") %
                     boost::to_lower_copy(EnumToString(shell.getL())) %
                     shell.getSize() % 1.0;
      for (const GaussianPrimitive& gaussian : shell) {
        outFile << boost::format("%22.10e %22.10e\n") % gaussian.decay() %
                       gaussian.contraction();
      }
    }
    outFile << " \n";
  }

}  // namespace xtp

void MoldenWriter::WriteFile(const std::string& filename,
                             const Orbitals& orbitals) {
  if (orbitals.hasDFTbasisName()) {
    _bs.Load(orbitals.getDFTbasisName());
    _basis.Fill(_bs, orbitals.QMAtoms());
  } else {
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

}  // namespace xtp
}  // namespace votca