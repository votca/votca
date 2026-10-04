/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE aopotential_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aopotential.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/polarsite.h"
#include "xtp_libint2.h"
#include <votca/tools/eigenio_matrixmarket.h>
using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aopotential_test)

BOOST_AUTO_TEST_CASE(aopotentials_test) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  AOMultipole esp;
  esp.FillPotential(aobasis, orbitals.QMAtoms());
  Eigen::MatrixXd esp_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/esp_ref.mm");
  bool check_esp = esp.Matrix().isApprox(esp_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_esp, 1);
  if (!check_esp) {
    std::cout << "esp Ref" << endl;
    std::cout << esp_ref << endl;
    std::cout << "esp" << endl;
    std::cout << esp.Matrix() << endl;
  }

  ECPBasisSet ecps;
  ecps.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/ecp.xml");
  ECPAOBasis ecpbasis;
  ecpbasis.Fill(ecps, orbitals.QMAtoms());
  AOECP ecp;
  ecp.FillPotential(aobasis, ecpbasis);
  Eigen::MatrixXd ecp_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/ecp_ref.mm");

  bool check_ecp = ecp.Matrix().isApprox(ecp_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_ecp, 1);
  if (!check_ecp) {
    std::cout << "ecp Ref" << endl;
    std::cout << ecp_ref << endl;
    std::cout << "ecp" << endl;
    std::cout << ecp.Matrix() << endl;
  }

  StaticSegment seg("", 0);
  seg.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/aopotential/polarsite.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites;
  for (const StaticSite& site : seg) {
    externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
  }
  AOMultipole dip;
  dip.FillPotential(aobasis, externalsites);

  Eigen::MatrixXd dip_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/dip_ref.mm");
  bool dip_check = dip_ref.isApprox(dip.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(dip_check, 1);
  if (!dip_check) {
    std::cout << "dip Ref" << endl;
    std::cout << dip_ref << endl;
    std::cout << "Dip" << endl;
    std::cout << dip.Matrix() << endl;
  }

  StaticSegment seg2("", 0);
  seg2.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                    "/aopotential/polarsite2.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites2;
  for (const StaticSite& site : seg2) {
    externalsites2.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
  }

  AOMultipole quad;
  quad.FillPotential(aobasis, externalsites2);

  Eigen::MatrixXd quad_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/quad_ref.mm");
  bool quad_check = quad_ref.isApprox(quad.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(quad_check, 1);
  if (!quad_check) {
    std::cout << "Quad Ref" << endl;
    std::cout << quad_ref << endl;
    std::cout << "Quad" << endl;
    std::cout << quad.Matrix() << endl;
  }

  AOPlanewave planewave;
  std::vector<Eigen::Vector3d> kpoints = {
      {1, 1, 1}, {2, 1, 1}, {-1, -1, -1}, {-2, -1, -1}};
  Eigen::MatrixXd planewave_ref = Eigen::MatrixXd::Zero(17, 17);
  planewave_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/planewave_ref.mm");

  planewave.FillPotential(aobasis, kpoints);
  Eigen::MatrixXd planewave_real = planewave.Matrix().real();
  bool planewave_check = planewave_ref.isApprox(planewave_real, 1e-4);

  BOOST_CHECK_EQUAL(planewave_check, 1);
  if (!planewave_check) {
    std::cout << "planewave Ref real" << endl;
    std::cout << planewave_ref << endl;
    std::cout << "planewave real" << endl;
    std::cout << planewave.Matrix().real() << endl;
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(aomultipole_comparison) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  {
    StaticSegment seg("", 0);
    seg.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                     "/aopotential/polarsite_ao.mps");

    std::vector<std::unique_ptr<StaticSite> > externalsites;
    for (const StaticSite& site : seg) {
      externalsites.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site)));
    }
    AOMultipole dip;
    dip.FillPotential(aobasis, externalsites);

    Eigen::MatrixXd dip_ref = Eigen::MatrixXd::Zero(17, 17);

    double a = 0.1;                            // this is in a0
    double mag_d = seg[0].getDipole().norm();  // this is in e * a0
    const Eigen::Vector3d dir_d = seg[0].getDipole().normalized();
    const Eigen::Vector3d A = seg[0].getPos() + 0.5 * a * dir_d;
    const Eigen::Vector3d B = seg[0].getPos() - 0.5 * a * dir_d;
    double qA = mag_d / a;
    double qB = -qA;
    StaticSite site1 = StaticSite(0, "", A);
    site1.setCharge(qA);
    StaticSite site2 = StaticSite(1, "", B);
    site2.setCharge(qB);
    std::vector<std::unique_ptr<StaticSite> > externalsites_mono;
    externalsites_mono.push_back(
        std::unique_ptr<StaticSite>(new StaticSite(site1)));
    externalsites_mono.push_back(
        std::unique_ptr<StaticSite>(new StaticSite(site2)));
    AOMultipole mono2;
    mono2.FillPotential(aobasis, externalsites_mono);

    bool dip_check = mono2.Matrix().isApprox(dip.Matrix(), 1e-4);
    BOOST_CHECK_EQUAL(dip_check, 1);
    if (!dip_check) {
      std::cout << "mono2 Ref" << endl;
      std::cout << mono2.Matrix() << endl;
      std::cout << "Dip" << endl;
      std::cout << dip.Matrix() << endl;
    }
  }

  {

    StaticSegment seg2("", 0);
    seg2.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                      "/aopotential/polarsite_ao2.mps");

    std::vector<std::unique_ptr<StaticSite> > externalsites2;
    for (const StaticSite& site : seg2) {
      externalsites2.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site)));
    }
    AOMultipole quad;
    quad.FillPotential(aobasis, externalsites2);

    std::vector<std::unique_ptr<StaticSite> > externalsites_mono6;
    const Eigen::Matrix3d components = seg2[0].CalculateCartesianMultipole();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(components);
    double a = 2 * 0.01;
    for (Index i = 0; i < 3; i++) {
      double q = es.eigenvalues()[i] / (a * a);
      const Eigen::Vector3d vec1 =
          seg2[0].getPos() + 0.5 * a * es.eigenvectors().col(i);
      const Eigen::Vector3d vec2 =
          seg2[0].getPos() - 0.5 * a * es.eigenvectors().col(i);
      StaticSite site1 = StaticSite(0, "", vec1);
      site1.setCharge(q);
      StaticSite site2 = StaticSite(1, "", vec2);
      site2.setCharge(q);
      externalsites_mono6.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site1)));
      externalsites_mono6.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site2)));
    }

    AOMultipole mono6;
    mono6.FillPotential(aobasis, externalsites_mono6);

    bool quad_check = mono6.Matrix().isApprox(quad.Matrix(), 1e-4);
    BOOST_CHECK_EQUAL(quad_check, 1);
    if (!quad_check) {
      std::cout << "mono6 Ref" << endl;
      std::cout << mono6.Matrix() << endl;
      std::cout << "Quad" << endl;
      std::cout << quad.Matrix() << endl;
      std::cout << "diff" << endl;
      std::cout << mono6.Matrix().cwiseQuotient(quad.Matrix()) << endl;
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(large_l_test) {
  libint2::initialize();
  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/C2.xyz");

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/G.xml");
  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);

  Index dftbasissize = 18;

  StaticSegment seg2("", 0);
  seg2.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                    "/aopotential/polarsite_all.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites2;
  for (const StaticSite& site : seg2) {
    externalsites2.push_back(std::make_unique<StaticSite>(site));
  }

  AOMultipole esp;
  esp.FillPotential(dftbasis, externalsites2);
  Eigen::MatrixXd esp_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/esp_ref_l.mm");
  bool check_esp = esp.Matrix().isApprox(esp_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_esp, 1);
  if (!check_esp) {
    cout << "esp ref" << endl;
    cout << esp_ref << endl;
    cout << "esp result" << endl;
    cout << esp.Matrix() << endl;
  }

  ECPBasisSet ecps;
  ecps.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/ecp.xml");
  ECPAOBasis ecpbasis;
  ecpbasis.Fill(ecps, mol);
  AOECP ecp;
  ecp.FillPotential(dftbasis, ecpbasis);

  Eigen::MatrixXd ecp_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/ecp_ref_l.mm");
  bool check_ecp = ecp.Matrix().isApprox(ecp_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_ecp, 1);
  if (!check_ecp) {
    cout << "ecp ref" << endl;
    cout << ecp_ref << endl;
    cout << "ecp result" << endl;
    cout << ecp.Matrix() << endl;
  }

  AOPlanewave planewave;
  std::vector<Eigen::Vector3d> kpoints = {
      {1, 1, 1}, {2, 1, 1}, {-1, -1, -1}, {-2, -1, -1}};
  Eigen::MatrixXd planewave_ref =
      Eigen::MatrixXd::Zero(dftbasissize, dftbasissize);
  planewave_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/planewave_ref_l.mm");

  planewave.FillPotential(dftbasis, kpoints);
  bool planewave_check =
      planewave_ref.isApprox(planewave.Matrix().real(), 1e-4);
  BOOST_CHECK_EQUAL(planewave_check, 1);
  if (!planewave_check) {
    std::cout << "planewave Ref real" << endl;
    std::cout << planewave_ref << endl;
    std::cout << "planewave real" << endl;
    std::cout << planewave.Matrix().real() << endl;
  }
  libint2::finalize();
}

// All sites are integrated in one pass over the shell pairs. Sites of
// different rank in one batch must give exactly the sum of their individual
// matrices, and the nuclear potential must equal the sum over atoms.
BOOST_AUTO_TEST_CASE(batched_sites_equal_sum_of_single_sites) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  std::srand(4);
  std::vector<std::unique_ptr<StaticSite> > sites;
  for (Index i = 0; i < 9; ++i) {
    auto site = std::make_unique<StaticSite>(
        i, "C", Eigen::Vector3d(4.0 * Eigen::Vector3d::Random()));
    Vector9d m = Vector9d::Zero();
    m(0) = 0.3 * Eigen::VectorXd::Random(1)(0);
    const Index rank = i % 3;  // charges, dipoles and quadrupoles mixed
    if (rank > 0) {
      m.segment<3>(1) = 0.2 * Eigen::Vector3d::Random();
    }
    if (rank > 1) {
      m.segment<5>(4) = 0.1 * Eigen::VectorXd::Random(5);
    }
    site->setMultipole(m, rank);
    sites.push_back(std::move(site));
  }
  AOMultipole all;
  all.FillPotential(aobasis, sites);

  Eigen::MatrixXd sum =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  for (const auto& site : sites) {
    std::vector<std::unique_ptr<StaticSite> > one;
    one.push_back(std::make_unique<StaticSite>(*site));
    AOMultipole single;
    single.FillPotential(aobasis, one);
    sum += single.Matrix();
  }
  BOOST_CHECK_SMALL((all.Matrix() - sum).cwiseAbs().maxCoeff(),
                    1e-13 * sum.cwiseAbs().maxCoeff());

  AOMultipole nuclei;
  nuclei.FillPotential(aobasis, orbitals.QMAtoms());
  Eigen::MatrixXd nuc_sum =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  for (const QMAtom& atom : orbitals.QMAtoms()) {
    std::vector<std::unique_ptr<StaticSite> > one;
    one.push_back(
        std::make_unique<StaticSite>(atom, double(atom.getNuccharge())));
    AOMultipole single;
    single.FillPotential(aobasis, one);
    nuc_sum += single.Matrix();
  }
  BOOST_CHECK_SMALL((nuclei.Matrix() - nuc_sum).cwiseAbs().maxCoeff(),
                    1e-13 * nuc_sum.cwiseAbs().maxCoeff());
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(permanent_plus_induced_equals_all) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  std::srand(7);
  std::vector<std::unique_ptr<StaticSite> > sites;
  for (Index i = 0; i < 12; ++i) {
    const Eigen::Vector3d pos = 4.0 * Eigen::Vector3d::Random();
    Vector9d m = Vector9d::Zero();
    m(0) = 0.3 * Eigen::VectorXd::Random(1)(0);
    const Index rank = i % 3;
    if (rank > 0) {
      m.segment<3>(1) = 0.2 * Eigen::Vector3d::Random();
    }
    if (rank > 1) {
      m.segment<5>(4) = 0.1 * Eigen::VectorXd::Random(5);
    }
    if (i % 4 == 3) {  // a static site among the polar ones
      auto site = std::make_unique<StaticSite>(i, "C", pos);
      site->setMultipole(m, rank);
      sites.push_back(std::move(site));
      continue;
    }
    auto site = std::make_unique<PolarSite>(i, "C", pos);
    site->setMultipole(m, rank);
    if (i % 4 != 2) {  // some polar sites without induced dipole
      site->setInduced_Dipole(0.1 * Eigen::Vector3d::Random());
    }
    sites.push_back(std::move(site));
  }
  AOMultipole all;
  all.FillPotential(aobasis, sites);
  AOMultipole permanent;
  permanent.FillPotential(aobasis, sites, AOMultipole::Moments::Permanent);
  AOMultipole induced;
  induced.FillPotential(aobasis, sites, AOMultipole::Moments::Induced);
  const Eigen::MatrixXd sum = permanent.Matrix() + induced.Matrix();
  BOOST_CHECK_SMALL((all.Matrix() - sum).cwiseAbs().maxCoeff(),
                    1e-13 * all.Matrix().cwiseAbs().maxCoeff());
  BOOST_CHECK_GT(induced.Matrix().cwiseAbs().maxCoeff(),
                 1e-3 * all.Matrix().cwiseAbs().maxCoeff());

  // without induced dipoles the permanent part is everything
  for (auto& site : sites) {
    if (auto* polar = dynamic_cast<PolarSite*>(site.get())) {
      polar->setInduced_Dipole(Eigen::Vector3d::Zero());
    }
  }
  AOMultipole all0;
  all0.FillPotential(aobasis, sites);
  BOOST_CHECK_SMALL((all0.Matrix() - permanent.Matrix()).cwiseAbs().maxCoeff(),
                    1e-13 * all0.Matrix().cwiseAbs().maxCoeff());
  libint2::finalize();
}

// Rank-1 sites use a recursion contracted with the dipole; a rank-2 site
// with zero quadrupole goes through the per-component dipole recursion. Both
// must agree, up to g functions on both sides of the integral.
BOOST_AUTO_TEST_CASE(contracted_dipole_equals_components) {
  libint2::initialize();
  auto Check = [](const AOBasis& aobasis) {
    std::srand(11);
    std::vector<std::unique_ptr<StaticSite> > dipoles;
    std::vector<std::unique_ptr<StaticSite> > as_rank2;
    for (Index i = 0; i < 6; ++i) {
      const Eigen::Vector3d pos = 3.0 * Eigen::Vector3d::Random();
      Vector9d m = Vector9d::Zero();
      m(0) = 0.3 * Eigen::VectorXd::Random(1)(0);
      m.segment<3>(1) = 0.5 * Eigen::Vector3d::Random();
      auto d = std::make_unique<StaticSite>(i, "C", pos);
      d->setMultipole(m, 1);
      dipoles.push_back(std::move(d));
      auto q = std::make_unique<StaticSite>(i, "C", pos);
      q->setMultipole(m, 2);
      as_rank2.push_back(std::move(q));
    }
    AOMultipole contracted;
    contracted.FillPotential(aobasis, dipoles);
    AOMultipole components;
    components.FillPotential(aobasis, as_rank2);
    const double scale = components.Matrix().cwiseAbs().maxCoeff();
    BOOST_CHECK_GT(scale, 1e-3);
    BOOST_CHECK_SMALL(
        (contracted.Matrix() - components.Matrix()).cwiseAbs().maxCoeff(),
        1e-13 * scale);
  };

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aopotential/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  Check(aobasis);

  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/C2.xyz");
  BasisSet gbasis;
  gbasis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aopotential/G.xml");
  AOBasis gaobasis;
  gaobasis.Fill(gbasis, mol);
  Check(gaobasis);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
