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

#define BOOST_TEST_MODULE chargeresponsefit_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/chargeresponsefit.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(chargeresponsefit_test)

namespace {

const std::string kBasis =
    std::string(XTP_TEST_DATA_FOLDER) + "/hirshfeldpartition/3-21G.xml";

QMMolecule Water() {
  QMMolecule mol("water", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/chargeresponse/water.xyz");
  return mol;
}

// xtpdft options; extra is inserted into <dftpackage> (e.g. a field)
votca::tools::Property DFTOptions(const std::string& extra = "") {
  const std::string file = "chargeresponsefit_test.xml";
  {
    std::ofstream xml(file);
    xml << "<dftpackage>\n<spin>1</spin>\n<name>xtp</name>\n"
        << "<charge>0</charge>\n"
        << "<functional>XC_HYB_GGA_XC_PBEH</functional>\n"
        << "<basisset>" << kBasis << "</basisset>\n"
        << "<auxbasisset>" << XTP_TEST_DATA_FOLDER
        << "/diabatization/aux-def2-svp.xml</auxbasisset>\n"
        << "<initial_guess>independent</initial_guess>\n"
        << extra << "<xtpdft>\n<screening_eps>1e-12</screening_eps>\n"
        << "<fock_matrix_reset>5</fock_matrix_reset>\n<convergence>\n"
        << "<energy>1e-12</energy>\n<method>DIIS</method>\n"
        << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
        << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
        << "<levelshift_end>0.2</levelshift_end>\n"
        << "<max_iterations>100</max_iterations>\n<error>1e-10</error>\n"
        << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
        << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
        << "<davidson_max_iter>50</davidson_max_iter>\n"
        << "</convergence>\n<integration_grid>fine</integration_grid>\n"
        << "<max_iterations>200</max_iterations>\n</xtpdft>\n"
        << "</dftpackage>\n";
  }
  votca::tools::Property prop;
  prop.LoadFromXML(file);
  return prop.get("dftpackage");
}

struct Partition {
  QMMolecule mol = Water();
  AOBasis basis;
  Vxc_Grid grid;
  std::vector<HirshfeldPartition::AtomicReference> atoms;
  Partition() {
    BasisSet bs;
    bs.Load(kBasis);
    basis.Fill(bs, mol);
    grid.GridSetup("fine", mol, basis);
    votca::tools::Property options = DFTOptions();
    Logger log;
    DFTEngine dft;
    dft.setLogger(&log);
    dft.Initialize(options);
    atoms = dft.HirshfeldReferences(mol);
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(weight_matrices_are_the_hirshfeld_weights) {
  libint2::initialize();
  const Partition p;
  const std::vector<Eigen::MatrixXd> W =
      ChargeResponseFit::WeightMatrices(p.atoms, p.basis, p.grid);
  BOOST_REQUIRE_EQUAL(W.size(), 3);
  AOOverlap overlap;
  overlap.Fill(p.basis);
  Eigen::MatrixXd sum = -overlap.Matrix();
  for (Index k = 0; k < 3; ++k) {
    const Eigen::MatrixXd ref =
        HirshfeldPartition::BuildWeightMatrix(p.atoms, k, p.basis, p.grid);
    BOOST_CHECK_SMALL((W[std::size_t(k)] - ref).cwiseAbs().maxCoeff(), 1e-12);
    sum += W[std::size_t(k)];
  }
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(atomic_polarizabilities_partition_the_dipole) {
  // sum_k alpha_k = -int r drho + sum_k R_k int w_k drho, for any density
  libint2::initialize();
  const Partition p;
  const std::vector<Eigen::MatrixXd> W =
      ChargeResponseFit::WeightMatrices(p.atoms, p.basis, p.grid);
  const Index n = p.basis.AOBasisSize();
  std::array<Eigen::MatrixXd, 3> dP;
  for (Eigen::MatrixXd& m : dP) {
    const Eigen::MatrixXd r = Eigen::MatrixXd::Random(n, n);
    m = 0.5 * (r + r.transpose());
  }
  std::vector<Eigen::Vector3d> positions;
  for (const QMAtom& atom : p.mol) {
    positions.push_back(atom.getPos());
  }
  const std::vector<Eigen::Matrix3d> alpha =
      ChargeResponseFit::AtomicPolarizabilities(p.atoms, positions, p.basis,
                                                p.grid, dP);
  AODipole dipole;
  dipole.setCenter(Eigen::Vector3d::Zero());
  dipole.Fill(p.basis);
  Eigen::Matrix3d sum = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d ref = Eigen::Matrix3d::Zero();
  for (Index c = 0; c < 3; ++c) {
    for (Index a = 0; a < 3; ++a) {
      ref(a, c) = -dipole.Matrix()[std::size_t(a)]
                       .cwiseProduct(dP[std::size_t(c)])
                       .sum();
    }
    for (Index k = 0; k < 3; ++k) {
      ref.col(c) += positions[std::size_t(k)] *
                    W[std::size_t(k)].cwiseProduct(dP[std::size_t(c)]).sum();
    }
  }
  for (const Eigen::Matrix3d& a : alpha) {
    sum += a;
  }
  BOOST_CHECK_SMALL((sum - ref).cwiseAbs().maxCoeff(),
                    1e-5 * ref.cwiseAbs().maxCoeff());
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(fit_of_water) {
  libint2::initialize();
  Logger log;
  log.setReportLevel(votca::Log::error);
  const votca::tools::Property options = DFTOptions();
  Orbitals orb;
  orb.QMAtoms() = Water();
  orb.setChargeAndSpin(0, 1);
  const ChargeResponseFit::Result r =
      ChargeResponseFit::Fit(options, orb, 1e-3, log);

  // finite differences reproduce the symmetry of the linear response
  BOOST_CHECK_SMALL(r.pi_asymmetry, 1e-5);
  BOOST_CHECK_SMALL(r.pi.rowwise().sum().cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK(r.pi_eigenvalues(1) > 0.0);
  BOOST_CHECK_SMALL((r.alpha_ct + r.alpha_res - r.alpha).norm(), 1e-12);
  BOOST_CHECK_CLOSE(r.res_weights.sum(), 1.0, 1e-10);
  // the two hydrogens are equivalent
  BOOST_CHECK_CLOSE(r.pi(1, 1), r.pi(2, 2), 1e-4);
  BOOST_CHECK_CLOSE(r.res_weights(1), r.res_weights(2), 1e-4);

  // independent check of the field responses: alpha_zz from SCFs in a
  // uniform field (xtpdft externalfield)
  const double F = 1e-3;
  std::array<double, 2> mu;
  for (std::size_t s = 0; s < 2; ++s) {
    const double f = (s == 0) ? F : -F;
    votca::tools::Property field = DFTOptions(
        "<externalfield>0 0 " + std::to_string(f) + "</externalfield>\n");
    DFTEngine dft;
    dft.setLogger(&log);
    dft.Initialize(field);
    Orbitals o;
    o.QMAtoms() = Water();
    o.setChargeAndSpin(0, 1);
    BOOST_REQUIRE(dft.Evaluate(o));
    BasisSet bs;
    bs.Load(kBasis);
    AOBasis basis;
    basis.Fill(bs, o.QMAtoms());
    AODipole dipole;
    dipole.setCenter(o.QMAtoms().getPos());
    dipole.Fill(basis);
    mu[s] =
        -dipole.Matrix()[2].cwiseProduct(o.DensityMatrixGroundState()).sum();
  }
  const double alpha_zz = std::abs(mu[0] - mu[1]) / (2 * F);
  BOOST_CHECK_CLOSE(r.alpha(2, 2), alpha_zz, 1e-2);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(inversion_symmetry) {
  // ethylene: the responses of half the carbons and hydrogens by symmetry
  // equal the full set; water has no inversion centre
  libint2::initialize();
  const std::string file = "chargeresponsefit_ethylene.xyz";
  {
    std::ofstream xyz(file);
    xyz << "6\nethylene\n"
        << "C 0.0 0.0 0.6695\nC 0.0 0.0 -0.6695\n"
        << "H 0.0 0.9289 1.2321\nH 0.0 -0.9289 1.2321\n"
        << "H 0.0 0.9289 -1.2321\nH 0.0 -0.9289 -1.2321\n";
  }
  QMMolecule mol("ethylene", 0);
  mol.LoadFromFile(file);
  BasisSet bs;
  bs.Load(kBasis);
  AOBasis basis;
  basis.Fill(bs, mol);
  const ChargeResponseFit::Inversion inv =
      ChargeResponseFit::FindInversion(mol, basis);
  BOOST_REQUIRE_EQUAL(inv.atom_map.size(), 6);
  BOOST_CHECK_EQUAL(inv.atom_map[0], 1);
  BOOST_CHECK_EQUAL(inv.atom_map[2], 5);
  AOOverlap overlap;
  overlap.Fill(basis);
  BOOST_CHECK_SMALL(
      (inv.T * overlap.Matrix() * inv.T.transpose() - overlap.Matrix())
          .cwiseAbs()
          .maxCoeff(),
      1e-12);
  {
    AOBasis water_basis;
    water_basis.Fill(bs, Water());
    BOOST_CHECK(ChargeResponseFit::FindInversion(Water(), water_basis)
                    .atom_map.empty());
  }

  Logger log;
  log.setReportLevel(votca::Log::error);
  const votca::tools::Property options = DFTOptions();
  std::array<ChargeResponseFit::Result, 2> r;
  for (std::size_t i = 0; i < 2; ++i) {
    Orbitals orb;
    orb.QMAtoms() = mol;
    orb.setChargeAndSpin(0, 1);
    r[i] = ChargeResponseFit::Fit(options, orb, 1e-3, log, i == 0);
  }
  BOOST_CHECK_SMALL((r[0].pi - r[1].pi).cwiseAbs().maxCoeff(),
                    1e-6 * r[1].pi.cwiseAbs().maxCoeff());
  BOOST_CHECK_SMALL((r[0].alpha_res - r[1].alpha_res).norm(),
                    1e-6 * r[1].alpha_res.norm());
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
