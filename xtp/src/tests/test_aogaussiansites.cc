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
#include "votca/xtp/dftengine.h"
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

BOOST_AUTO_TEST_CASE(quadrupole_point_limit) {
  libint2::initialize();
  Orbitals orb;
  AOBasis basis = Basis(orb);
  const Eigen::Vector3d R = orb.QMAtoms().getPos() + Eigen::Vector3d(3, 2, -2);
  StaticSite site(0, "H", R);
  Vector9d Q = Vector9d::Zero();
  Q.segment<5>(4) << 0.7, -0.3, 0.4, 0.25, -0.6;
  site.setMultipole(Q, 2);
  std::vector<std::unique_ptr<StaticSite>> sites;
  sites.push_back(std::make_unique<StaticSite>(site));
  AOMultipole m;
  m.FillPotential(basis, sites);
  AOGaussianSites::Source src{R, AOGaussianSites::kPointWidth, 0.0,
                              Eigen::Vector3d::Zero(),
                              site.CalculateCartesianMultipole()};
  const Eigen::MatrixXd g = AOGaussianSites::Potential(basis, {src});
  // O(width^2) from the point: 1e-5 for the width of the gradients
  BOOST_CHECK_SMALL((-g - m.Matrix()).norm() / m.Matrix().norm(), 2e-5);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(gradient_finite_difference) {
  libint2::initialize();
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    libint2::finalize();
    return;
  }
  Orbitals orb;
  AOBasis basis = Basis(orb);
  const Eigen::Vector3d c = orb.QMAtoms().getPos();
  Eigen::Matrix3d theta;
  theta << 0.3, 0.1, -0.2, 0.1, -0.5, 0.15, -0.2, 0.15, 0.2;
  std::vector<AOGaussianSites::Source> sources = {
      {c + Eigen::Vector3d(2.5, 1.0, -1.5), AOGaussianSites::kPointWidth, 0.6,
       Eigen::Vector3d(0.2, -0.4, 0.3), theta},
      {c + Eigen::Vector3d(-1.0, 2.0, 1.0), 1.1, -0.4,
       Eigen::Vector3d(-0.3, 0.1, 0.5)}};
  const Index n = basis.AOBasisSize();
  std::srand(5);
  Eigen::MatrixXd D = Eigen::MatrixXd::Random(n, n);
  D = 0.5 * (D + D.transpose()).eval();
  Eigen::MatrixXd site_grad;
  const Eigen::MatrixXd grad =
      AOGaussianSites::Gradient(basis, sources, D, site_grad);

  const double h = 1e-4;
  BasisSet bs;
  bs.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  auto energy = [&](const QMMolecule& mol,
                    const std::vector<AOGaussianSites::Source>& src) {
    AOBasis b;
    b.Fill(bs, mol);
    return AOGaussianSites::Potential(b, src).cwiseProduct(D).sum();
  };
  double max_err = 0.0;
  double scale = 0.0;
  for (Index a = 0; a < orb.QMAtoms().size(); ++a) {
    for (Index x = 0; x < 3; ++x) {
      QMMolecule plus = orb.QMAtoms();
      QMMolecule minus = orb.QMAtoms();
      Eigen::Vector3d shift = Eigen::Vector3d::Zero();
      shift(x) = h;
      plus[a].Translate(shift);
      minus[a].Translate(-shift);
      const double fd =
          (energy(plus, sources) - energy(minus, sources)) / (2 * h);
      max_err = std::max(max_err, std::abs(fd - grad(a, x)));
      scale = std::max(scale, std::abs(fd));
    }
  }
  for (std::size_t k = 0; k < sources.size(); ++k) {
    for (Index x = 0; x < 3; ++x) {
      std::vector<AOGaussianSites::Source> plus = sources;
      std::vector<AOGaussianSites::Source> minus = sources;
      plus[k].pos(x) += h;
      minus[k].pos(x) -= h;
      const double fd =
          (energy(orb.QMAtoms(), plus) - energy(orb.QMAtoms(), minus)) /
          (2 * h);
      max_err = std::max(max_err, std::abs(fd - site_grad(Index(k), x)));
      scale = std::max(scale, std::abs(fd));
    }
  }
  // translational invariance: atoms and sources together
  const double net =
      (grad.colwise().sum() + site_grad.colwise().sum()).cwiseAbs().maxCoeff();
  std::cout << "Gaussian-site gradient: max |analytic - FD| " << max_err
            << " (largest component " << scale << "), net " << net << std::endl;
  BOOST_CHECK_SMALL(max_err / scale, 1e-5);
  BOOST_CHECK_SMALL(net / scale, 1e-8);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(point_gradient_finite_difference) {
  libint2::initialize();
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    libint2::finalize();
    return;
  }
  // gradient of Tr(D V) with V the exact point-multipole potential of
  // AOMultipole (charge, dipole, quadrupole), for a site close to the
  // molecule and one further out
  Orbitals orb;
  AOBasis basis = Basis(orb);
  const Eigen::Vector3d c = orb.QMAtoms().getPos();
  std::vector<StaticSite> sites;
  {
    StaticSite near(
        0, "H", orb.QMAtoms()[1].getPos() + Eigen::Vector3d(1.4, 0.9, -0.6));
    Vector9d Q;
    Q << 0.4, 0.2, -0.3, 0.25, 0.5, -0.2, 0.3, 0.15, -0.4;
    near.setMultipole(Q, 2);
    sites.push_back(near);
    StaticSite far(1, "H", c + Eigen::Vector3d(-3.0, 2.5, 1.5));
    Q << -0.6, -0.1, 0.4, 0.2, -0.3, 0.1, 0.2, -0.25, 0.35;
    far.setMultipole(Q, 2);
    sites.push_back(far);
  }
  auto sources_of = [](const std::vector<StaticSite>& s) {
    std::vector<AOGaussianSites::Source> src;
    for (const StaticSite& site : s) {
      src.push_back({site.getPos(), 0.0, site.getCharge(), site.getDipole(),
                     site.CalculateCartesianMultipole()});
    }
    return src;
  };
  const Index n = basis.AOBasisSize();
  std::srand(9);
  Eigen::MatrixXd D = Eigen::MatrixXd::Random(n, n);
  D = 0.5 * (D + D.transpose()).eval();
  Eigen::MatrixXd site_grad;
  const Eigen::MatrixXd grad =
      AOGaussianSites::PointGradient(basis, sources_of(sites), D, site_grad);

  BasisSet bs;
  bs.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  // Tr(D Potential) = -Tr(D V_AOMultipole)
  auto energy = [&](const QMMolecule& mol, const std::vector<StaticSite>& s) {
    AOBasis b;
    b.Fill(bs, mol);
    std::vector<std::unique_ptr<StaticSite>> ptrs;
    for (const StaticSite& site : s) {
      ptrs.push_back(std::make_unique<StaticSite>(site));
    }
    AOMultipole m;
    m.FillPotential(b, ptrs);
    return -m.Matrix().cwiseProduct(D).sum();
  };
  const double h = 1e-4;
  double max_err = 0.0;
  double scale = 0.0;
  for (Index a = 0; a < orb.QMAtoms().size(); ++a) {
    for (Index x = 0; x < 3; ++x) {
      QMMolecule plus = orb.QMAtoms();
      QMMolecule minus = orb.QMAtoms();
      Eigen::Vector3d shift = Eigen::Vector3d::Zero();
      shift(x) = h;
      plus[a].Translate(shift);
      minus[a].Translate(-shift);
      const double fd = (energy(plus, sites) - energy(minus, sites)) / (2 * h);
      max_err = std::max(max_err, std::abs(fd - grad(a, x)));
      scale = std::max(scale, std::abs(fd));
    }
  }
  for (std::size_t k = 0; k < sites.size(); ++k) {
    for (Index x = 0; x < 3; ++x) {
      std::vector<StaticSite> plus = sites;
      std::vector<StaticSite> minus = sites;
      Eigen::Vector3d shift = Eigen::Vector3d::Zero();
      shift(x) = h;
      plus[k].Translate(shift);
      minus[k].Translate(-shift);
      const double fd =
          (energy(orb.QMAtoms(), plus) - energy(orb.QMAtoms(), minus)) /
          (2 * h);
      max_err = std::max(max_err, std::abs(fd - site_grad(Index(k), x)));
      scale = std::max(scale, std::abs(fd));
    }
  }
  std::cout << "point-multipole gradient: max |analytic - FD| " << max_err
            << " (largest component " << scale << ")" << std::endl;
  BOOST_CHECK_SMALL(max_err / scale, 1e-6);
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
