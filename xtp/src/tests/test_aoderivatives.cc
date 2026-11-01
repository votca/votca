/*
 * Copyright 2009-2024 The VOTCA Development Team (http://www.votca.org)
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

// Finite-difference tests of the nuclear derivatives of AO integrals in
// libint2_derivative_calls.cc for H2 with 3-21G: overlap, kinetic, nuclear
// attraction, two-centre Coulomb metric and three-centre RI integrals, each
// also checked for translational invariance. Skipped if libint2 lacks
// derivative integrals for the operator.

#include "xtp_libint2.h"
#include <stdexcept>

// Full libint2 API (Engine, Shell, BraKet) for the deriv_order=0 reference
// of three_center_derivative_finite_difference via ComputeAO3cBlock;
// xtp_libint2.h only provides initialize()/finalize().
// <libint2/statics_definition.h> must not be included: it defines libint2's
// static tables, which libvotca_xtp (libint2_calls.cc) already provides.
#define LIBINT2_CONSTEXPR_STATICS 0
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wcpp"
#endif
#include <libint2.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aoderivatives_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

// Defined in libint2_derivative_calls.cc (no public header).
namespace votca {
namespace xtp {
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis& aobasis, const QMMolecule& mol);
using ThreeCenterDerivative = std::array<std::vector<Eigen::MatrixXd>, 3>;
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis& auxbasis, const AOBasis& dftbasis);
// Declared (not in any header) in libint2_calls.cc; used here only to
// build the deriv_order=0 finite-difference reference.
std::vector<Eigen::MatrixXd> ComputeAO3cBlock(const libint2::Shell& auxshell,
                                              const AOBasis& dftbasis,
                                              libint2::Engine& engine);
}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_SUITE(aoderivatives_test)

// H2 along z: two atoms exercise the derivatives w.r.t. both shell centres,
// including the s1 != s2 transposed blocks.
QMMolecule BuildH2(double bond_length_angstrom) {
  QMMolecule mol(" ", 0);
  std::string xyz_content =
      "2\n\n"
      "H 0.0 0.0 0.0\n"
      "H 0.0 0.0 " +
      std::to_string(bond_length_angstrom) + "\n";
  std::string tmp_path = "/tmp/xtp_test_h2_tmp.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();
  mol.LoadFromFile(tmp_path);
  return mol;
}

BOOST_AUTO_TEST_CASE(overlap_derivative_finite_difference) {
  libint2::initialize();

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // finite-difference step, Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis aobasis0;
  aobasis0.Fill(basisset, mol0);

  std::vector<AOMatrixDerivative> deriv;
  try {
    deriv = ComputeOverlapDerivatives(aobasis0);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING overlap_derivative_finite_difference: " << e.what()
              << std::endl;
    libint2::finalize();
    return;
  }
  // deriv[a][xyz]: derivative w.r.t. atom a. For H2 along z only z is
  // nonzero, and translational invariance gives deriv[0][2] = -deriv[1][2].
  BOOST_CHECK_SMALL((deriv[0][2] + deriv[1][2]).cwiseAbs().maxCoeff(), 1e-8);

  // Central finite difference: displace atom 1 (index 1, the second H)
  // along z by +-h, recompute overlap, compare to analytic deriv[1][2].
  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basisset, mol_plus);
  AOOverlap overlap_plus;
  overlap_plus.Fill(aobasis_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basisset, mol_minus);
  AOOverlap overlap_minus;
  overlap_minus.Fill(aobasis_minus);

  Eigen::MatrixXd finite_diff_deriv =
      (overlap_plus.Matrix() - overlap_minus.Matrix()) / (2.0 * h);
  // h is in Angstrom, the analytic derivative per Bohr.
  constexpr double kAngstromPerBohr = 0.52917721090380;
  finite_diff_deriv *= kAngstromPerBohr;

  // observed relative deviation ~2.4e-5, the O(h^2) error at h = 1e-4
  bool matches = finite_diff_deriv.isApprox(deriv[1][2], 1e-4);
  if (!matches) {
    std::cout << "Analytic dS/dz(atom1):\n" << deriv[1][2] << std::endl;
    std::cout << "Finite-difference dS/dz(atom1):\n"
              << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails with a large (not small numerical) "
                 "discrepancy, check the buffer-ordering assumption flagged "
                 "in libint2_derivative_calls.cc first."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// As overlap_derivative_finite_difference, for the kinetic energy
// integrals (same code path, different libint2 operator).
BOOST_AUTO_TEST_CASE(kinetic_derivative_finite_difference) {
  libint2::initialize();

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // finite-difference step, Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis aobasis0;
  aobasis0.Fill(basisset, mol0);

  std::vector<AOMatrixDerivative> deriv;
  try {
    deriv = ComputeKineticDerivatives(aobasis0);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING kinetic_derivative_finite_difference: " << e.what()
              << std::endl;
    libint2::finalize();
    return;
  }

  // Translational invariance: the kinetic integral depends only on the
  // separation of the two centres.
  BOOST_CHECK_SMALL((deriv[0][2] + deriv[1][2]).cwiseAbs().maxCoeff(), 1e-8);

  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basisset, mol_plus);
  AOKinetic kinetic_plus;
  kinetic_plus.Fill(aobasis_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basisset, mol_minus);
  AOKinetic kinetic_minus;
  kinetic_minus.Fill(aobasis_minus);

  Eigen::MatrixXd finite_diff_deriv =
      (kinetic_plus.Matrix() - kinetic_minus.Matrix()) / (2.0 * h);
  // per Bohr, see overlap test
  constexpr double kAngstromPerBohr = 0.52917721090380;
  finite_diff_deriv *= kAngstromPerBohr;

  // same tolerance as the overlap test
  bool matches = finite_diff_deriv.isApprox(deriv[1][2], 1e-4);
  if (!matches) {
    std::cout << "Analytic dT/dz(atom1):\n" << deriv[1][2] << std::endl;
    std::cout << "Finite-difference dT/dz(atom1):\n"
              << finite_diff_deriv << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// Two-centre Coulomb metric derivative d(P|Q)/dR, using the 3-21G basis as
// the "auxiliary" basis.
BOOST_AUTO_TEST_CASE(coulomb_metric_derivative_finite_difference) {
  libint2::initialize();

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // finite-difference step, Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis aobasis0;
  aobasis0.Fill(basisset, mol0);

  std::vector<AOMatrixDerivative> deriv;
  try {
    deriv = ComputeCoulombMetricDerivatives(aobasis0);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING coulomb_metric_derivative_finite_difference: "
              << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  // Translational invariance: (P|Q) depends only on the separation.
  BOOST_CHECK_SMALL((deriv[0][2] + deriv[1][2]).cwiseAbs().maxCoeff(), 1e-8);

  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basisset, mol_plus);
  AOCoulomb coulomb_plus;
  coulomb_plus.Fill(aobasis_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basisset, mol_minus);
  AOCoulomb coulomb_minus;
  coulomb_minus.Fill(aobasis_minus);

  Eigen::MatrixXd finite_diff_deriv =
      (coulomb_plus.Matrix() - coulomb_minus.Matrix()) / (2.0 * h);
  constexpr double kAngstromPerBohr = 0.52917721090380;
  finite_diff_deriv *= kAngstromPerBohr;

  bool matches = finite_diff_deriv.isApprox(deriv[1][2], 1e-4);
  if (!matches) {
    std::cout << "Analytic d(P|Q)/dz(atom1):\n" << deriv[1][2] << std::endl;
    std::cout << "Finite-difference d(P|Q)/dz(atom1):\n"
              << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails with a large (not small numerical) "
                 "discrepancy, the buffer-count hypothesis for this "
                 "operator/braket combination (see libint2_derivative_"
                 "calls.cc) is the first thing to check."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// Three-centre RI integral derivative d(mu,nu|P)/dR, with H2/3-21G as both
// DFT and auxiliary basis (only the integral derivative is tested). A large
// discrepancy would point at the buffer-to-centre order (aux, col, row).
BOOST_AUTO_TEST_CASE(three_center_derivative_finite_difference) {
  libint2::initialize();

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // finite-difference step, Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis dftbasis0;
  dftbasis0.Fill(basisset, mol0);
  AOBasis auxbasis0;
  auxbasis0.Fill(basisset, mol0);  // same basis set used for both

  std::vector<ThreeCenterDerivative> deriv;
  try {
    deriv = ComputeThreeCenterDerivatives(auxbasis0, dftbasis0);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING three_center_derivative_finite_difference: "
              << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  // Translational invariance: with only two atoms the three-centre sum rule
  // reduces to d/dR_0 + d/dR_1 = 0 for every aux function and xyz.
  Index n_aux_bf = auxbasis0.AOBasisSize();
  for (Index k = 0; k < n_aux_bf; ++k) {
    BOOST_CHECK_SMALL((deriv[0][2][k] + deriv[1][2][k]).cwiseAbs().maxCoeff(),
                      1e-8);
  }

  // Finite-difference reference from (mu,nu|P) via ComputeAO3cBlock
  // (deriv_order = 0) at the displaced geometries.
  auto build_3c_tensor = [&](const AOBasis& auxbasis, const AOBasis& dftbasis) {
    std::vector<Eigen::MatrixXd> tensor(
        auxbasis.AOBasisSize(),
        Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));
    std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
    std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();
    libint2::Engine engine(
        libint2::Operator::coulomb,
        std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
        static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
    engine.set(libint2::BraKet::xs_xx);
    for (Index aux = 0; aux < auxbasis.getNumofShells(); ++aux) {
      std::vector<Eigen::MatrixXd> block =
          ComputeAO3cBlock(auxshells[aux], dftbasis, engine);
      Index aux_start = auxshell2bf[aux];
      for (size_t i = 0; i < block.size(); ++i) {
        tensor[aux_start + static_cast<Index>(i)] = block[i];
      }
    }
    return tensor;
  };

  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis dftbasis_plus;
  dftbasis_plus.Fill(basisset, mol_plus);
  AOBasis auxbasis_plus;
  auxbasis_plus.Fill(basisset, mol_plus);
  std::vector<Eigen::MatrixXd> tensor_plus =
      build_3c_tensor(auxbasis_plus, dftbasis_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis dftbasis_minus;
  dftbasis_minus.Fill(basisset, mol_minus);
  AOBasis auxbasis_minus;
  auxbasis_minus.Fill(basisset, mol_minus);
  std::vector<Eigen::MatrixXd> tensor_minus =
      build_3c_tensor(auxbasis_minus, dftbasis_minus);

  constexpr double kAngstromPerBohr = 0.52917721090380;
  bool all_match = true;
  for (Index k = 0; k < n_aux_bf; ++k) {
    Eigen::MatrixXd finite_diff_deriv =
        (tensor_plus[k] - tensor_minus[k]) / (2.0 * h) * kAngstromPerBohr;
    bool matches = finite_diff_deriv.isApprox(deriv[1][2][k], 1e-4);
    if (!matches) {
      all_match = false;
      std::cout << "Mismatch for aux function " << k << ":\n";
      std::cout << "Analytic d(mu,nu|P)/dz(atom1):\n"
                << deriv[1][2][k] << std::endl;
      std::cout << "Finite-difference:\n" << finite_diff_deriv << std::endl;
    }
  }
  if (!all_match) {
    std::cout << "NOTE: if this fails with a large (not small numerical) "
                 "discrepancy, check the buffer-count/ordering hypothesis "
                 "flagged in libint2_derivative_calls.cc first."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(all_match, true);

  libint2::finalize();
}

// Nuclear attraction integral derivative d(V_ne)/dR, a three-centre
// derivative per point charge (two basis centres and the nucleus).
BOOST_AUTO_TEST_CASE(nuclear_attraction_derivative_finite_difference) {
  libint2::initialize();

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis aobasis0;
  aobasis0.Fill(basisset, mol0);

  // Some libint2 builds lack derivatives for the nuclear attraction
  // operator only; ComputeNuclearAttractionDerivatives then throws, and the
  // test is skipped.
  std::vector<AOMatrixDerivative> deriv;
  try {
    deriv = ComputeNuclearAttractionDerivatives(aobasis0, mol0);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING nuclear_attraction_derivative_finite_difference: "
              << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  // Translational invariance: V_ne depends only on relative positions of
  // basis centres and nuclei.
  BOOST_CHECK_SMALL((deriv[0][2] + deriv[1][2]).cwiseAbs().maxCoeff(), 1e-8);

  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basisset, mol_plus);
  AOMultipole esp_plus;
  esp_plus.FillPotential(aobasis_plus, mol_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basisset, mol_minus);
  AOMultipole esp_minus;
  esp_minus.FillPotential(aobasis_minus, mol_minus);

  Eigen::MatrixXd finite_diff_deriv =
      (esp_plus.Matrix() - esp_minus.Matrix()) / (2.0 * h);
  constexpr double kAngstromPerBohr = 0.52917721090380;
  finite_diff_deriv *= kAngstromPerBohr;

  bool matches = finite_diff_deriv.isApprox(deriv[1][2], 1e-4);
  if (!matches) {
    std::cout << "Analytic d(V_ne)/dz(atom1):\n" << deriv[1][2] << std::endl;
    std::cout << "Finite-difference d(V_ne)/dz(atom1):\n"
              << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails with a large (not small numerical) "
                 "discrepancy, check the buffer-count/ordering hypothesis "
                 "for this operator (9 buffers, [shell1 atom][shell2 atom]"
                 "[point-charge atom]) flagged in "
                 "libint2_derivative_calls.cc first. If the SIGN looks "
                 "flipped specifically, check the explicit negation "
                 "against AOMultipole::FillPotential's own convention "
                 "(also documented there)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
