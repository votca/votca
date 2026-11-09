/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#include <cmath>
#include <fstream>
#include <stdexcept>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/diabatwriter.h"

namespace votca {
namespace xtp {

void DiabatWriter::WriteFile(const std::string& filename,
                             const Orbitals& orbitals, QMStateType type,
                             Index nstates) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    throw std::runtime_error("DiabatWriter: could not open " + filename);
  }
  Write(out, orbitals, type, nstates);
}

void DiabatWriter::Write(std::ostream& out, const Orbitals& orbitals,
                         QMStateType type, Index nstates) const {

  if (orbitals.hasUnrestrictedOrbitals() || orbitals.getSpin() != 1) {
    throw std::runtime_error(
        "DiabatWriter: only restricted closed-shell BSE is supported.");
  }
  if (type != QMStateType::Singlet && type != QMStateType::Triplet) {
    throw std::runtime_error(
        "DiabatWriter: only singlet or triplet BSE states can be written.");
  }
  const bool singlet = (type == QMStateType::Singlet);
  if ((singlet && !orbitals.hasBSESinglets()) ||
      (!singlet && !orbitals.hasBSETriplets())) {
    throw std::runtime_error("DiabatWriter: orbitals contain no BSE " +
                             type.ToLongString() + " states.");
  }
  const tools::EigenSystem& es =
      singlet ? orbitals.BSESinglets() : orbitals.BSETriplets();

  const bool tda = orbitals.getTDAApprox();
  const Index nbasis = orbitals.getBasisSetSize();
  // doubly occupied levels (valence only if ECPs are used, consistent with
  // the fchk written from the same orbitals)
  const Index nocc = orbitals.getLumo();
  const Index vmin = orbitals.getBSEvmin();
  const Index cmin = orbitals.getBSEcmin();
  const Index vtotal = orbitals.getBSEvmax() - vmin + 1;
  const Index ctotal = orbitals.getBSEcmax() - cmin + 1;
  if (orbitals.getBSEvmax() != nocc - 1 || cmin != nocc) {
    throw std::runtime_error(
        "DiabatWriter: BSE window does not start at HOMO/LUMO.");
  }
  if (es.eigenvectors().rows() != vtotal * ctotal) {
    throw std::runtime_error(
        "DiabatWriter: BSE eigenvector size does not match the BSE window.");
  }
  if (!tda && es.eigenvectors2().rows() != vtotal * ctotal) {
    throw std::runtime_error(
        "DiabatWriter: full BSE requested but no anti-resonant (Y) part "
        "stored.");
  }

  Index nwrite = es.eigenvalues().size();
  if (nstates > 0 && nstates < nwrite) {
    nwrite = nstates;
  }

  Eigen::VectorXd osc = Eigen::VectorXd::Zero(nwrite);
  if (singlet && orbitals.hasTransitionDipoles()) {
    Eigen::VectorXd all = orbitals.Oscillatorstrengths();
    const Index n = std::min<Index>(nwrite, all.size());
    osc.head(n) = all.head(n);
  }

  const double e0 = orbitals.getDFTTotalEnergy();
  const double spin_factor = 1.0 / std::sqrt(2.0);
  const std::string rule(60, '-');

  out << FormatIdentifier << "\n";
  out << "# Written by VOTCA-XTP (orb2diabat) from BSE "
      << (singlet ? "singlets" : "triplets") << ", "
      << (tda ? "Tamm-Dancoff" : "full BSE") << ", window MO " << vmin + 1
      << ".." << orbitals.getBSEcmax() + 1 << "\n";
  out << "Method: " << (tda ? "TDA" : "TDDFT") << "\n";
  out << "Orbitals are not spin-polarized.\n";
  out << "Number of electrons: " << 2 * nocc << "\n";
  out << "Number of alpha electrons: " << nocc << "\n";
  out << "Number of beta electrons: " << nocc << "\n";
  out << "Number of occupied orbitals (alpha): " << nocc << "\n";
  out << "Number of virtual orbitals (alpha): " << nbasis - nocc << "\n";
  out << "Number of occupied orbitals (beta): " << nocc << "\n";
  out << "Number of virtual orbitals (beta): " << nbasis - nocc << "\n";
  out << boost::format(
             "Total energy of the reference ground state (a.u.): "
             "%.12E\n") %
             e0;
  out << "\nState Information:\n\n";

  for (Index i = 0; i < nwrite; ++i) {
    const double omega = es.eigenvalues()(i);
    Eigen::VectorXd x = spin_factor * es.eigenvectors().col(i);
    Eigen::VectorXd y;
    if (!tda) {
      y = spin_factor * es.eigenvectors2().col(i);
    }
    const double nx = x.squaredNorm();
    const double ny = tda ? 0.0 : y.squaredNorm();

    out << "\n  State index: " << i + 1 << "\n\n";
    out << boost::format("  Total energy (a.u.): %.12E\n") % (e0 + omega);
    out << boost::format("  Excitation energy (eV): %16.12f\n") %
               (omega * tools::conv::hrt2ev);
    out << "  Multiplicity: " << (singlet ? "Singlet" : "Triplet") << "\n";
    out << boost::format("  <S^2>: %9.6f\n") % (singlet ? 0.0 : 2.0);
    out << boost::format("  Strength: %12.6f\n") % osc(i);
    out << boost::format("  Normalization factor: %9.6f\n") % (2.0 * (nx - ny));
    out << "  Amplitudes:\n  " << rule << "\n";
    // VOTCA compound index of a (v,c) pair is v * ctotal + c
    for (Index v = 0; v < vtotal; ++v) {
      for (Index c = 0; c < ctotal; ++c) {
        out << boost::format("  %6d   -->  %5d   % .15E  alpha\n") %
                   (vmin + v + 1) % (cmin + c + 1) % x(v * ctotal + c);
      }
    }
    if (!tda) {
      for (Index v = 0; v < vtotal; ++v) {
        for (Index c = 0; c < ctotal; ++c) {
          out << boost::format("  %6d   <--  %5d   % .15E  alpha\n") %
                     (vmin + v + 1) % (cmin + c + 1) % y(v * ctotal + c);
        }
      }
    }
    out << "  " << rule << "\n";
    out << boost::format("    X: %12.8f   Y: %11.8f\n") % (2.0 * nx) %
               (2.0 * ny);
  }
  out << "\n";
}

}  // namespace xtp
}  // namespace votca
