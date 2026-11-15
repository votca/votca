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
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/excitonfollow.h"
#include "votca/xtp/qmmolecule.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd MixedAOOverlap(const Orbitals& a, const Orbitals& b) {
  const QMMolecule& ma = a.QMAtoms();
  const QMMolecule& mb = b.QMAtoms();
  if (ma.size() != mb.size()) {
    throw std::runtime_error("MixedAOOverlap: different numbers of atoms");
  }
  for (Index i = 0; i < ma.size(); ++i) {
    if (ma[i].getElement() != mb[i].getElement()) {
      throw std::runtime_error("MixedAOOverlap: different atoms");
    }
  }
  if (a.getDFTbasisName() != b.getDFTbasisName()) {
    throw std::runtime_error("MixedAOOverlap: different basis sets");
  }
  if (a.getBasisSetSize() != b.getBasisSetSize()) {
    throw std::runtime_error("MixedAOOverlap: basis size mismatch");
  }
  return ComputeMixedOverlap(a.getDftBasis(), b.getDftBasis());
}

namespace {
// AO representation C_v M^T C_c^T of the exciton column (X or Y)
Eigen::MatrixXd ExcitonAO(const Orbitals& orb, const Eigen::VectorXd& col) {
  const Index vmin = orb.getBSEvmin();
  const Index nv = orb.getBSEvmax() - vmin + 1;
  const Index cmin = orb.getBSEcmin();
  const Index nc = orb.getBSEcmax() - cmin + 1;
  const Eigen::MatrixXd& C = orb.MOs().eigenvectors();
  Eigen::Map<const Eigen::MatrixXd> m(col.data(), nc, nv);
  return C.middleCols(vmin, nv) * m.transpose() *
         C.middleCols(cmin, nc).transpose();
}

const tools::EigenSystem& Excitons(const Orbitals& orb, QMStateType type) {
  if (type == QMStateType::Singlet) {
    if (!orb.hasBSESinglets()) {
      throw std::runtime_error("ExcitonOverlaps: no BSE singlets");
    }
    return orb.BSESinglets();
  } else if (type == QMStateType::Triplet) {
    if (!orb.hasBSETriplets()) {
      throw std::runtime_error("ExcitonOverlaps: no BSE triplets");
    }
    return orb.BSETriplets();
  }
  throw std::runtime_error("ExcitonOverlaps: singlets or triplets only");
}
}  // namespace

Eigen::VectorXd ExcitonOverlaps(const Orbitals& ref, const QMState& ref_state,
                                const Orbitals& cur, QMStateType type) {
  const tools::EigenSystem& eref = Excitons(ref, ref_state.Type());
  const tools::EigenSystem& ecur = Excitons(cur, type);
  const Index k = ref_state.StateIdx();
  if (k < 0 || k >= eref.eigenvectors().cols()) {
    throw std::runtime_error("ExcitonOverlaps: reference state not stored");
  }
  const Eigen::MatrixXd S = MixedAOOverlap(ref, cur);  // ref x cur
  // reference exciton carried to the current basis, then to the current
  // occupied x virtual MO pairs: T = C_v^T S^T A_ref S C_c (nv x nc)
  const Index vmin = cur.getBSEvmin();
  const Index nv = cur.getBSEvmax() - vmin + 1;
  const Index cmin = cur.getBSEcmin();
  const Index nc = cur.getBSEcmax() - cmin + 1;
  const Eigen::MatrixXd& C = cur.MOs().eigenvectors();
  const Eigen::MatrixXd SCv = S * C.middleCols(vmin, nv);
  const Eigen::MatrixXd SCc = S * C.middleCols(cmin, nc);
  auto project = [&](const Eigen::VectorXd& col) {
    const Eigen::MatrixXd A = ExcitonAO(ref, col);
    // stored like the exciton columns: (c, v), column major
    const Eigen::MatrixXd T = (SCv.transpose() * A * SCc).transpose();
    return Eigen::VectorXd(
        Eigen::Map<const Eigen::VectorXd>(T.data(), T.size()));
  };
  Eigen::VectorXd ov =
      ecur.eigenvectors().transpose() * project(eref.eigenvectors().col(k));
  const bool full = !ref.getTDAApprox() && !cur.getTDAApprox() &&
                    eref.eigenvectors2().cols() > 0 &&
                    ecur.eigenvectors2().cols() > 0;
  if (full) {
    ov -=
        ecur.eigenvectors2().transpose() * project(eref.eigenvectors2().col(k));
  }
  return ov.cwiseAbs();
}

}  // namespace xtp
}  // namespace votca
