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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE environmentscreening_cr_test

#include <boost/test/unit_test.hpp>

#include "votca/xtp/chargeresponsesolver.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/qmmolecule.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(environmentscreening_cr_test)

namespace {

QMMolecule Methane() {
  QMMolecule mol("methane", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  return mol;
}

AOBasis MakeBasis(const std::string& name, const QMMolecule& mol) {
  BasisSet bs;
  bs.Load(name);
  AOBasis basis;
  basis.Fill(bs, mol);
  return basis;
}

// a water-like three-site molecule at offset, as a charge-response segment
// (Pi of a bond-incidence form) or a Thole segment
PolarSegment Molecule(Index id, const Eigen::Vector3d& offset, bool cr,
                      double width) {
  PolarSegment seg("mol", id);
  const std::vector<Eigen::Vector3d> pos = {Eigen::Vector3d(0, 0, 0.22),
                                            Eigen::Vector3d(0, 1.43, -0.89),
                                            Eigen::Vector3d(0, -1.43, -0.89)};
  const std::vector<std::string> el = {"O", "H", "H"};
  for (Index i = 0; i < 3; ++i) {
    PolarSite site(i, el[std::size_t(i)], pos[std::size_t(i)] + offset);
    site.setpolarization((i == 0 ? 5.5 : 2.0) * Eigen::Matrix3d::Identity());
    if (cr) {
      site.setChargeResponse(width, false);
    }
    seg.push_back(site);
  }
  if (cr) {
    Eigen::MatrixXd pi(3, 3);
    pi << 2.4, -1.2, -1.2, -1.2, 1.5, -0.3, -1.2, -0.3, 1.5;
    seg.setChargeResponse(ChargeResponse(el, pi));
  }
  return seg;
}

std::vector<Eigen::Vector3d> Positions(const std::vector<PolarSegment>& segs) {
  std::vector<Eigen::Vector3d> pos;
  for (const PolarSegment& seg : segs) {
    for (const PolarSite& site : seg) {
      pos.push_back(site.getPos());
    }
  }
  return pos;
}
}  // namespace

BOOST_AUTO_TEST_CASE(potential_and_field_agree) {
  // F (AuxFieldAtPoints) is minus the gradient of AuxPotentialAtPoints, for
  // Gaussian and point sites
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  const Eigen::Vector3d P = mol.getPos() + Eigen::Vector3d(3.5, -2.0, 1.5);
  for (double w : {0.0, 1.3}) {
    const Eigen::MatrixXd F = EnvironmentScreening::AuxFieldAtPoints(
        aux, {P}, std::vector<double>{w});
    const double h = 1e-3;
    for (Index k = 0; k < 3; ++k) {
      auto phi = [&](double t) {
        Eigen::Vector3d Q = P;
        Q(k) += t;
        return EnvironmentScreening::AuxPotentialAtPoints(aux, {Q}, {w});
      };
      const Eigen::MatrixXd grad =
          (-phi(2 * h) + 8 * phi(h) - 8 * phi(-h) + phi(-2 * h)) / (12 * h);
      BOOST_CHECK((F.col(k) + grad).cwiseAbs().maxCoeff() <
                  1e-8 * F.col(k).cwiseAbs().maxCoeff());
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(thole_limit) {
  // without charge-response segments the general kernel is the Thole one
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  std::vector<PolarSegment> segs = {
      Molecule(0, mol.getPos() + Eigen::Vector3d(6, 0, 0), false, 0.0),
      Molecule(1, mol.getPos() + Eigen::Vector3d(0, 6.5, 1), false, 0.0)};
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, Positions(segs));
  const Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(F.rows(), 0);
  const Eigen::MatrixXd B1 =
      EnvironmentScreening::ReactionFieldKernel(F, segs, 0.39);
  const Eigen::MatrixXd B2 =
      EnvironmentScreening::ChargeResponseKernel(F, Phi, segs, 0.39);
  BOOST_CHECK((B1 - B2).cwiseAbs().maxCoeff() <
              1e-10 * B1.cwiseAbs().maxCoeff());
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(kernel_is_the_region_response) {
  // c^T B c for an auxiliary density is twice the polarization energy of
  // the polar region (ChargeResponseSolver, iterative) in its potential and
  // field
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  const double w = 1.3;
  std::vector<PolarSegment> segs = {
      Molecule(0, mol.getPos() + Eigen::Vector3d(5, 0, 0), true, w),
      Molecule(1, mol.getPos() + Eigen::Vector3d(0, 5.5, 1), false, 0.0),
      Molecule(2, mol.getPos() + Eigen::Vector3d(-1, -1, 5.5), true, w)};
  const std::vector<double> widths =
      EnvironmentScreening::SiteWidths(segs, 0.0);
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, Positions(segs), widths);
  std::vector<Eigen::Vector3d> cr_pos;
  std::vector<double> cr_w;
  for (const PolarSegment& seg : segs) {
    for (const PolarSite& site : seg) {
      if (site.isChargeResponse()) {
        cr_pos.push_back(site.getPos());
        cr_w.push_back(site.CRWidth());
      }
    }
  }
  BOOST_REQUIRE_EQUAL(cr_pos.size(), 6);
  BOOST_CHECK_EQUAL(widths[0], w);
  const Eigen::MatrixXd Phi =
      EnvironmentScreening::AuxPotentialAtPoints(aux, cr_pos, cr_w);
  const Eigen::MatrixXd B =
      EnvironmentScreening::ChargeResponseKernel(F, Phi, segs, 0.39);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(aux.AOBasisSize());
  c(0) = 1.0;
  c(3) = -0.4;
  c(aux.AOBasisSize() - 1) = 0.7;
  // the same density as sources for the iterative solver
  const Eigen::VectorXd phi = Phi.transpose() * c;
  const Eigen::VectorXd E = F.transpose() * c;
  Index q = 0, s = 0;
  for (PolarSegment& seg : segs) {
    for (PolarSite& site : seg) {
      site.Reset();
      if (site.isChargeResponse()) {
        site.Phi_noE() = phi(q++);
      }
      site.V_noE() = -E.segment<3>(3 * s++);
    }
  }
  const eeInteractor ee(0.39);
  ChargeResponseSolver solver(segs, ee);
  BOOST_REQUIRE(solver.Solve(1e-13, 200, false).converged);
  // E_pol = 1/2 (dq . phi - mu . E)
  double e_pol = 0.0;
  q = 0;
  s = 0;
  for (const PolarSegment& seg : segs) {
    for (const PolarSite& site : seg) {
      if (site.isChargeResponse()) {
        e_pol += 0.5 * site.Induced_Charge() * phi(q++);
      }
      e_pol -= 0.5 * site.Induced_Dipole().dot(E.segment<3>(3 * s++));
    }
  }
  BOOST_CHECK_CLOSE(0.5 * c.dot(B * c), e_pol, 1e-8);
  BOOST_CHECK(e_pol < 0.0);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
