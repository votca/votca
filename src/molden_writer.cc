#include "votca/xtp/molden_writer.h"

#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

void Molden_Writer::writeAtoms(const Orbitals& orbitals,
                               std::ofstream& outFile) {
  for (auto& atom : orbitals.QMAtoms()) {
    Eigen::Vector3d pos = atom.getPos();
    outFile << boost::format("%-4s %5d %5d %22.12e %22.10e %22.10e\n") %
                   atom.getElement() % (atom.getId() + 1) %
                   atom.getPureNucCharge() % pos[0] % pos[1] % pos[2];
  }
}

void Molden_Writer::writeMOs(const Orbitals& orbitals, std::ofstream& outFile) {

  Eigen::VectorXd energies = orbitals.MOs().eigenvalues();

  OrbReorder reorder(_transpositions, _multipliers);

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

void Molden_Writer::writeBasisSet(const Orbitals& orbitals,
                                  std::ofstream& outFile) {
  if (orbitals.hasDFTbasisName()) {

    for (auto& atom : orbitals.QMAtoms()) {
      const Element& element = _bs.getElement(atom.getElement());
      // The 0 in the format string of the next line is meaningless it
      // is included for backwards compatibility of molden files
      outFile << boost::format("%4d 0 \n") % (atom.getId() + 1);
      for (const Shell& shell : element) {
        for (const char& subtype : shell.getType()) {
          // The 1.0 at the end of the next line is meaningless it
          // is included for backwards compatibility of molden files
          outFile << boost::format("%-3s %4d %3.1f \n") %
                         std::tolower(subtype, std::locale()) %
                         shell.getSize() % 1.0;
          for (const GaussianPrimitive& gaussian : shell) {
            outFile << boost::format("%22.10e %22.10e\n") % gaussian.decay() %
                           gaussian.Contractions()[FindLmax(
                               std::string(1, subtype))];
          }
        }
      }
      outFile << " \n";
    }

  } else {
    throw std::runtime_error(".orb file does not contain a basisset name");
  }
}

void Molden_Writer::WriteFile(const std::string& filename,
                              const Orbitals& orbitals) const {

  _bs.Load(orbitals.getDFTbasisName());
  _basis.Fill(_bs, orbitals.QMAtoms());

  std::ofstream outFile(filename);

  if (outFile.is_open()) {
    XTP_LOG(Log::info, _log) << "Writing data to " << filename << std::flush;
    // print Header
    outFile << "[Molden Format]\n";
    outFile << "[Title]\n";
    outFile << "Molden file created by VOTCA-XTP for basename: " << _job_name
            << "\n";
    outFile << " \n";

    outFile << "[Atoms] AU\n";
    writeAtoms(orbitals, outFile);

    outFile << "[GTO] \n";
    writeBasisSet(orbitals, outFile);

    // indicate spherical D F and G functions
    outFile << "[5D] \n[7F] \n[9G] \n";

    outFile << "[MO]\n";
    writeMOs(orbitals, outFile);

    XTP_LOG(Log::error, _log) << "Done parsing \n" << std::flush;
    return true;
  }
  return false;
}

}  // namespace xtp
}  // namespace votca