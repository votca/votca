/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/basisset.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sstream>
#include <votca/tools/eigenio_matrixmarket.h>

namespace votca {
namespace xtp {

/*
 * This function converts VOTCA's L enum to the gaussian equivalent.
 * Gaussian uses a minus sign to indicate spherical shells, but -1 for sp.
 */
Index GaussianWriter::toGaussianL(L l) const {
  switch (l) {
    case L::S:
      return 0;
    case L::P:
      return 1;
    default:
      return -1 * static_cast<Index>(l);
  }
}

std::string GaussianWriter::reorderedMOCoefficients(
    const Orbitals& orbitals) const {
  OrbReorder reorder(gaussianOrder, gaussianMultipliers, true);
  Eigen::MatrixXd moCoefficients = orbitals.MOs().eigenvectors();
  reorder.reorderOrbitals(moCoefficients, orbitals.SetupDftBasis());

  // put the reordered mos in a string
  std::stringstream mos_string;

  int temp_int = 1;
  for (Index i = 0; i < moCoefficients.rows(); ++i) {
    for (Index j = 0; j < moCoefficients.cols(); ++j) {
      mos_string << boost::format("%16.8e") % moCoefficients(j, i);
      if (temp_int % 5 == 0) {
        mos_string << "\n";
      }
      temp_int++;
    }
  }
  mos_string << ((temp_int - 1) % 5 == 0 ? "" : "\n");

  return mos_string.str();
}

std::string GaussianWriter::densityMatrixToString(const Orbitals& orbitals,
                                                  const QMState& state,
                                                  bool diff2gs) const {
  OrbReorder reorder(gaussianOrder, gaussianMultipliers, true);

  Eigen::MatrixXd density;
  if (state.Type().isExciton() && diff2gs) {
    std::array<Eigen::MatrixXd, 2> DMAT =
        orbitals.DensityMatrixExcitedState(state);
    density = DMAT[1] - DMAT[0];
  } else if ((state.Type().isKSState() || state.Type().isPQPState()) &&
             diff2gs) {
    density = orbitals.DensityMatrixKSstate(state);
  } else {
    density = orbitals.DensityMatrixFull(state);
  }

  reorder.reorderOperator(density, orbitals.SetupDftBasis());

  // put the reordered mos in a string
  std::stringstream density_string;

  int temp_int = 1;

  for (Index i = 0; i < density.rows(); ++i) {
    for (Index j = 0; j <= i; ++j) {
      density_string << boost::format("%16.8e") % density(i, j);
      if (temp_int % 5 == 0) {
        density_string << "\n";
      }
      temp_int++;
    }
  }
  density_string << ((temp_int - 1) % 5 == 0 ? "" : "\n");

  return density_string.str();
}

void GaussianWriter::WriteFile(const std::string& basename,
                               const Orbitals& orbitals, const QMState state,
                               bool diff2gs) const {
  if (!orbitals.hasDFTbasisName()) {
    throw std::runtime_error(".orb file does not contain a basisset name");
  }

  AOBasis basis = orbitals.SetupDftBasis();

  std::ofstream outFile(basename + ".fchk");

  if (outFile.is_open()) {
    XTP_LOG(Log::error, _log)
        << "Start writing to " << (basename + ".fchk") << std::flush;
    int temp_int;
    // job description
    outFile << basename << ", fchk created by VOTCA-XTP\n";
    outFile << "SP    RHF    " << orbitals.getDFTbasisName() << "\n";

    // clang-format off
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of atoms" % "I" % orbitals.QMAtoms().size();
    outFile << boost::format("%-43s%-2s%15d\n") % "Charge" % "I" % 0;
    outFile << boost::format("%-43s%-2s%15d\n") % "Multiplicity" % "I" % 1;
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of electrons" % "I" % (2*orbitals.getNumberOfAlphaElectrons());
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of alpha electrons" % "I" % orbitals.getNumberOfAlphaElectrons();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of beta electrons" % "I" % orbitals.getNumberOfAlphaElectrons();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of basis functions" % "I" % basis.AOBasisSize();
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of independent functions" % "I" % basis.AOBasisSize();
    // clang-format on
    // ATOMIC NUMBERS
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Atomic numbers" % "I" %
                   orbitals.QMAtoms().size();
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      outFile << boost::format("%12d") % atom.getElementNumber();
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // NUCLEAR CHARGES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Nuclear charges" % "R" %
                   orbitals.QMAtoms().size();
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      outFile << boost::format("%16.8e") % (double)atom.getNuccharge();
      if (temp_int % 5 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // CURRENT CARTESIAN COORDINATES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Current cartesian coordinates" % "R" %
                   (3 * orbitals.QMAtoms().size());
    temp_int = 1;
    for (const auto& atom : orbitals.QMAtoms()) {
      for (int i = 0; i < 3; ++i) {
        outFile << boost::format("%16.8e") % atom.getPos()(i);
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // NUMBER OF PRIMITIVE SHELLS
    outFile << boost::format("%-43s%-2s%15d\n") % "Number of primitive shells" %
                   "I" % basis.getNumberOfPrimitives();
    // NUMBER OF CONTRACTED SHELLS
    outFile << boost::format("%-43s%-2s%15d\n") %
                   "Number of contracted shells" % "I" % basis.getNumofShells();
    // PURE/CARTESIAN D
    outFile << boost::format("%-43s%-2s%15d\n") % "Pure/Cartesian d shells " %
                   "I" % 0;
    // PURE/CARTESIAN F
    outFile << boost::format("%-43s%-2s%15d\n") % "Pure/Cartesian f shells " %
                   "I" % 0;
    // HIGHEST ANGULAR MOMENTUM
    outFile << boost::format("%-43s%-2s%15d\n") % "Highest angular momentum " %
                   "I" % basis.getMaxL();
    // HIGHEST ANGULAR MOMENTUM
    outFile << boost::format("%-43s%-2s%15d\n") %
                   "Largest degree of contraction " % "I" % basis.getMaxNprim();
    // SHELL TYPES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Shell types" % "I" %
                   (basis.getNumofShells());
    temp_int = 1;
    for (const auto& shell : basis) {
      outFile << boost::format("%12d") % toGaussianL(shell.getL());
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // NR OF PRIMITIVES PER SHELL
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Number of primitives per shell" % "I" %
                   (basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      outFile << boost::format("%12d") % shell.getSize();
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // SHELL TO ATOM MAP
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Shell to atom map" %
                   "I" % (basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      // Gaussian indices start at 1, hence the + 1
      outFile << boost::format("%12d") % (shell.getAtomIndex() + 1);
      if (temp_int % 6 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 6 == 0 ? "" : "\n");
    // PRIMITIVE EXPONENTS
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Primitive exponents" %
                   "R" % basis.getNumberOfPrimitives();
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (const AOGaussianPrimitive& prim : shell) {
        outFile << boost::format("%16.8e") % prim.getDecay();
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // CONTRACTION COEFFICIENTS
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Contraction coefficients" % "R" %
                   basis.getNumberOfPrimitives();
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (const AOGaussianPrimitive& prim : shell) {
        outFile << boost::format("%16.8e") % prim.getContraction();
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // SHELL COORDINATES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Coordinates of each shell" % "R" %
                   (3 * basis.getNumofShells());
    temp_int = 1;
    for (const AOShell& shell : basis) {
      for (int i = 0; i < 3; ++i) {
        outFile << boost::format("%16.8e") % shell.getPos()(i);
        if (temp_int % 5 == 0) {
          outFile << "\n";
        }
        temp_int++;
      }
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // TOTAL ENERGY
    outFile << boost::format("%-43s%-2s%22.15e\n") % "Total Energy" % "R" %
                   orbitals.getDFTTotalEnergy();
    // ALPHA ORBITAL ENERGIES
    outFile << boost::format("%-43s%-2s N=  %10d\n") %
                   "Alpha Orbital Energies" % "R" %
                   orbitals.MOs().eigenvalues().size();
    temp_int = 1;
    for (Index i = 0; i < orbitals.MOs().eigenvalues().size(); ++i) {
      outFile << boost::format("%16.8e") % orbitals.MOs().eigenvalues()[i];
      if (temp_int % 5 == 0) {
        outFile << "\n";
      }
      temp_int++;
    }
    outFile << ((temp_int - 1) % 5 == 0 ? "" : "\n");
    // ALPHA ORBITAL ENERGIES
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Alpha MO coefficients" %
                   "R" %
                   (orbitals.MOs().eigenvalues().size() *
                    orbitals.MOs().eigenvalues().size());
    outFile << reorderedMOCoefficients(orbitals);
    // DENSITY MATRIX
    outFile << boost::format("%-43s%-2s N=  %10d\n") % "Total SCF Density" %
                   "R" %
                   ((orbitals.MOs().eigenvalues().size() *
                     (orbitals.MOs().eigenvalues().size() - 1)) /
                        2 +
                    orbitals.MOs().eigenvalues().size());
    outFile << densityMatrixToString(orbitals, state, diff2gs);

    XTP_LOG(Log::error, _log) << "Done writing \n" << std::flush;
  }
}

}  // namespace xtp
}  // namespace votca