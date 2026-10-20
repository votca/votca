/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aogaussiansites_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aogaussiansites.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/polarsite.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using namespace votca;

namespace {
AOBasis Basis(Orbitals& orbitals) {
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  return aobasis;
}

Eigen::MatrixXd PointSites(const AOBasis& basis, const Eigen::Vector3d& pos,
                           double q, const Eigen::Vector3d& dipole) {
  StaticSite site(0, "H", pos);
  Vector9d Q = Vector9d::Zero();
  Q(0) = q;
  Q.segment<3>(1) = dipole;
  site.setMultipole(Q, dipole.norm() > 0 ? 1 : 0);
  std::vector<std::unique_ptr<StaticSite>> sites;
  sites.push_back(std::make_unique<StaticSite>(site));
  AOMultipole m;
  m.FillPotential(basis, sites);
  return m.Matrix();
}
}  // namespace

BOOST_AUTO_TEST_SUITE(aogaussiansites_test)

BOOST_AUTO_TEST_CASE(point_limit) {
  libint2::initialize();
  Orbitals orb;
  AOBasis basis = Basis(orb);
  const Eigen::Vector3d R = orb.QMAtoms().getPos() + Eigen::Vector3d(3, 2, -2);
  // a very narrow Gaussian is a point: the electrons' Hamiltonian
  // (-Potential) equals AOMultipole's
  const double s = 1e-3;
  Eigen::MatrixXd q =
      AOGaussianSites::Potential(basis, {{R, s, 0.8, Eigen::Vector3d::Zero()}});
  Eigen::MatrixXd q_ref = PointSites(basis, R, 0.8, Eigen::Vector3d::Zero());
  BOOST_CHECK((-q).isApprox(q_ref, 1e-8));
  const Eigen::Vector3d mu(0.3, -0.5, 0.2);
  Eigen::MatrixXd d = AOGaussianSites::Potential(basis, {{R, s, 0.0, mu}});
  Eigen::MatrixXd d_ref = PointSites(basis, R, 0.0, mu);
  BOOST_CHECK((-d).isApprox(d_ref, 1e-7));
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(dipole_is_derivative_of_charge) {
  libint2::initialize();
  Orbitals orb;
  AOBasis basis = Basis(orb);
  const Eigen::Vector3d R = orb.QMAtoms().getPos() + Eigen::Vector3d(1, 1, 0.5);
  const double s = 1.3;  // overlaps the molecule
  const double h = 1e-3;
  for (Index c = 0; c < 3; ++c) {
    const Eigen::Vector3d e = Eigen::Vector3d::Unit(c);
    auto charge_at = [&](double t) {
      return AOGaussianSites::Potential(
          basis, {{R + t * e, s, 1.0, Eigen::Vector3d::Zero()}});
    };
    // 4-point central difference, error O(h^4)
    Eigen::MatrixXd fd = (-charge_at(2 * h) + 8 * charge_at(h) -
                          8 * charge_at(-h) + charge_at(-2 * h)) /
                         (12 * h);
    Eigen::MatrixXd d = AOGaussianSites::Potential(basis, {{R, s, 0.0, e}});
    BOOST_CHECK((d - fd).cwiseAbs().maxCoeff() <
                1e-8 * d.cwiseAbs().maxCoeff());
  }
  // the potential and gradient of a density are the traces with the same
  // matrices
  Eigen::MatrixXd D =
      Eigen::MatrixXd::Random(basis.AOBasisSize(), basis.AOBasisSize());
  D = 0.5 * (D + D.transpose()).eval();
  Eigen::MatrixXd pg =
      AOGaussianSites::PotentialAndGradientAt(basis, D, {R}, {s});
  Eigen::MatrixXd q =
      AOGaussianSites::Potential(basis, {{R, s, 1.0, Eigen::Vector3d::Zero()}});
  BOOST_CHECK_CLOSE(pg(0, 0), q.cwiseProduct(D).sum(), 1e-9);
  for (Index c = 0; c < 3; ++c) {
    Eigen::MatrixXd d = AOGaussianSites::Potential(
        basis, {{R, s, 0.0, Eigen::Vector3d::Unit(c)}});
    BOOST_CHECK_CLOSE(pg(0, 1 + c), d.cwiseProduct(D).sum(), 1e-9);
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(overlap_of_gaussians) {
  // one s primitive density chi^2 against a Gaussian site: closed form
  // erf(r/S)/r with S^2 = s^2 + 1/(2a) (both are normalized Gaussians)
  libint2::initialize();
  QMMolecule mol("H", 0);
  mol.push_back(QMAtom(0, "H", Eigen::Vector3d::Zero()));
  BasisSet bs;
  std::ofstream f("one_s.xml");
  f << "<basis name=\"one_s\"><element name=\"H\"><shell type=\"S\" "
       "scale=\"1.0\"><constant decay=\"0.8\"><contractions type=\"S\" "
       "factor=\"1.0\"/></constant></shell></element></basis>\n";
  f.close();
  bs.Load("one_s.xml");
  AOBasis basis;
  basis.Fill(bs, mol);
  BOOST_REQUIRE_EQUAL(basis.AOBasisSize(), 1);
  const double a = 0.8, s = 1.1;
  const Eigen::Vector3d R(0.7, -0.4, 1.2);
  Eigen::MatrixXd D = Eigen::MatrixXd::Ones(1, 1);
  Eigen::MatrixXd pg =
      AOGaussianSites::PotentialAndGradientAt(basis, D, {R}, {s});
  const double S = std::sqrt(s * s + 1.0 / (2.0 * a));
  const double r = R.norm();
  BOOST_CHECK_CLOSE(pg(0, 0), std::erf(r / S) / r, 1e-8);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
