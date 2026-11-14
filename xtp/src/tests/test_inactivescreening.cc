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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE inactivescreening_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/inactivescreening.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using votca::Index;
namespace tools = votca::tools;

BOOST_AUTO_TEST_SUITE(inactivescreening_test)

namespace {

// Water with one hydrogen active (PbE through the O-H bond), as in
// test_activedftengine, and the embedded MOs of that run.
Orbitals EmbeddedWater() {
  Orbitals orb;
  orb.ReadFromCpt(std::string(XTP_TEST_DATA_FOLDER) +
                  "/activedftengine/molecule.orb");
  std::ofstream xml("dftengine.xml");
  xml << "<dftpackage>\n"
      << "<spin>1</spin>\n<name>xtp</name>\n<charge>0</charge>\n"
      << "<basisset>" << XTP_TEST_DATA_FOLDER
      << "/diabatization/def2-svp.xml</basisset>\n"
      << "<initial_guess>atom</initial_guess>\n"
      << "<functional>XC_HYB_GGA_XC_PBEH</functional>\n"
      << "<xtpdft>\n<screening_eps>1e-9</screening_eps>\n"
      << "<fock_matrix_reset>5</fock_matrix_reset>\n"
      << "<convergence>\n<energy>1e-9</energy>\n<method>DIIS</method>\n"
      << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
      << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
      << "<levelshift_end>0.2</levelshift_end>\n"
      << "<max_iterations>100</max_iterations>\n<error>1e-9</error>\n"
      << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
      << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
      << "<davidson_max_iter>50</davidson_max_iter>\n</convergence>\n"
      << "<integration_grid>xcoarse</integration_grid>\n"
      << "<max_iterations>100</max_iterations>\n"
      << "<dft_in_dft>\n<activeatoms>1</activeatoms>\n"
      << "<threshold>0.4</threshold>\n<levelshift>10000</levelshift>\n"
      << "<truncate_basis>false</truncate_basis>\n"
      << "<truncation_threshold>1e-4</truncation_threshold>\n"
      << "</dft_in_dft>\n</xtpdft>\n</dftpackage>\n";
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");
  Logger log;
  DFTEngine dft;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  orb.setChargeAndSpin(0, 1);
  BOOST_REQUIRE(dft.EvaluateActiveRegion(orb));
  return orb;
}

AOBasis AuxBasis(const Orbitals& orb) {
  BasisSet bs;
  bs.Load(std::string(XTP_TEST_DATA_FOLDER) +
          "/diabatization/aux-def2-svp.xml");
  AOBasis aux;
  aux.Fill(bs, orb.QMAtoms());
  return aux;
}

double MaxDiff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
  return (a - b).cwiseAbs().maxCoeff();
}

}  // namespace

BOOST_AUTO_TEST_CASE(static_screening_of_the_inactive_region) {
  libint2::initialize();
  Orbitals orb = EmbeddedWater();
  const AOBasis& dft = orb.getDftBasis();
  const AOBasis aux = AuxBasis(orb);

  const Index nocc = orb.getNumberOfAlphaElectrons();
  const Index nvirt = orb.MOs().eigenvectors().cols() - nocc;
  const Index nA = orb.getNumOfActiveElectrons() / 2;
  const tools::EigenSystem& emb = orb.getEmbeddedMOs();

  // inactive space: count, idempotency, orthogonality to the active MOs
  // (up to the leakage of the finite level shift, mu = 1e4 Ha here)
  const InactiveSpace B = InactiveScreening::BuildSpace(orb);
  BOOST_CHECK_EQUAL(B.size(), nocc - nA);
  BOOST_CHECK_SMALL(B.idempotency_error, 1e-10);
  AOOverlap S;
  S.Fill(dft);
  BOOST_CHECK_SMALL(
      (emb.eigenvectors().leftCols(nA).transpose() * S.Matrix() * B.occupied)
          .cwiseAbs()
          .maxCoeff(),
      1e-4);

  // reference: static response of the full calculation, and of the
  // embedded active orbitals (the inactive levels sit at the top of the
  // embedded MOs, shifted by the projector, and are left out)
  const Eigen::MatrixXd& C = orb.MOs().eigenvectors();
  const Eigen::VectorXd& e = orb.MOs().eigenvalues();
  const Eigen::MatrixXd X_full =
      InactiveScreening::Response(aux, dft, C.leftCols(nocc), e.head(nocc),
                                  C.rightCols(nvirt), e.tail(nvirt));
  const Eigen::MatrixXd X_A = InactiveScreening::Response(
      aux, dft, emb.eigenvectors().leftCols(nA), emb.eigenvalues().head(nA),
      emb.eigenvectors().middleCols(nA, nvirt),
      emb.eigenvalues().segment(nA, nvirt));
  const Index naux = X_full.rows();
  const Eigen::MatrixXd one = Eigen::MatrixXd::Identity(naux, naux);
  const Eigen::MatrixXd W_full = InverseSPD(one + X_full);
  // static W with the active integrals dressed by R
  auto dressed_W = [&](const Eigen::MatrixXd& R) {
    const Eigen::MatrixXd D = EnvironmentScreening::DressingMatrix(R);
    return Eigen::MatrixXd(D * InverseSPD(one + D * X_A * D) * D);
  };
  const double w_none = MaxDiff(InverseSPD(one + X_A), W_full);

  for (auto mode : {InactiveScreeningInput::Mode::Exact,
                    InactiveScreeningInput::Mode::Block}) {
    const InactiveScreeningInput in = InactiveScreening::Prepare(orb, mode);
    BOOST_CHECK_EQUAL(in.n_inactive, nocc - nA);
    const Eigen::MatrixXd Y =
        InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9);
    // one occupied orbital per chunk gives the same Y
    BOOST_CHECK_SMALL(
        MaxDiff(Y, InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 1.0)),
        1e-12);
    const Eigen::MatrixXd R = InactiveScreening::ReactionField(Y);
    const Eigen::VectorXd r = SymmetricEigen(R, false).values;
    BOOST_CHECK_LT(r.maxCoeff(), 1e-7);
    BOOST_CHECK_GT(r.minCoeff(), -1.0);
    const double w = MaxDiff(dressed_W(R), W_full);
    BOOST_TEST_MESSAGE(InactiveScreeningInput::Name(mode)
                       << ": max |W - W_full| = " << w
                       << " (without inactive screening " << w_none << ")");
    if (mode == InactiveScreeningInput::Mode::Exact) {
      // the full static W by construction
      BOOST_CHECK_SMALL(w, 1e-10);
    } else {
      // misses the active-inactive Fock coupling, which is large for a cut
      // through an O-H bond (0.65 Ha): 0.0081 vs 0.148 without
      BOOST_CHECK_LT(w, 0.1 * w_none);
    }
    // a zero environment changes nothing
    BOOST_CHECK_SMALL(MaxDiff(R, InactiveScreening::ReactionField(
                                     Y, Eigen::MatrixXd::Zero(naux, naux))),
                      1e-12);
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(kernels_from_the_stored_integrals) {
  // Kernels (one integral pass, X_A from the integrals GW holds) gives the
  // Kernel of each frequency (separate passes, X_A rebuilt)
  libint2::initialize();
  Orbitals orb = EmbeddedWater();
  const AOBasis& dft = orb.getDftBasis();
  const AOBasis aux = AuxBasis(orb);
  const Index nocc = orb.getNumberOfAlphaElectrons();
  const Index nvirt = orb.MOs().eigenvectors().cols() - nocc;
  const Index nA = orb.getNumOfActiveElectrons() / 2;
  const tools::EigenSystem& emb = orb.getEmbeddedMOs();
  const Index homo = nA - 1;
  const std::vector<double> w{1e-3, 0.5};
  // RPA range up to the full virtuals, and up to all embedded levels (the
  // inactive levels shifted by the projector sit on top and are left out)
  for (Index rpamax : {homo + nvirt, Index(emb.eigenvalues().size()) - 1}) {
    TCMatrix_gwbse Mmn;
    Mmn.Initialize(aux.AOBasisSize(), 0, rpamax, 0, rpamax);
    Mmn.Fill(aux, dft, emb.eigenvectors());
    const Eigen::VectorXd energies = emb.eigenvalues().head(rpamax + 1);
    for (auto mode : {InactiveScreeningInput::Mode::Exact,
                      InactiveScreeningInput::Mode::Block}) {
      const InactiveScreeningInput in = InactiveScreening::Prepare(orb, mode);
      const std::vector<Eigen::MatrixXd> Y = InactiveScreening::Kernels(
          aux, dft, in, Mmn, energies, homo, 0, rpamax, 4e9, w);
      BOOST_REQUIRE_EQUAL(Y.size(), w.size());
      for (std::size_t f = 0; f < w.size(); ++f) {
        const Eigen::MatrixXd ref =
            InactiveScreening::Kernel(aux, dft, in, emb, homo, 4e9, w[f]);
        BOOST_CHECK_SMALL(MaxDiff(Y[f], ref) / ref.cwiseAbs().maxCoeff(),
                          1e-10);
      }
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(inactive_response_in_the_rpa) {
  // exact mode at the plasmon-pole fitting points: the RPA of the embedded
  // integrals plus the inactive response reproduces the full epsilon, also
  // after the auxiliary index was rotated (as the PPM and the BSE do)
  libint2::initialize();
  Orbitals orb = EmbeddedWater();
  const AOBasis& dft = orb.getDftBasis();
  const AOBasis aux = AuxBasis(orb);
  const Index nocc = orb.getNumberOfAlphaElectrons();
  const Index nvirt = orb.MOs().eigenvectors().cols() - nocc;
  const Index nA = orb.getNumOfActiveElectrons() / 2;
  const tools::EigenSystem& emb = orb.getEmbeddedMOs();
  Logger log;

  TCMatrix_gwbse Mfull;
  Mfull.Initialize(aux.AOBasisSize(), 0, nocc + nvirt - 1, 0, nocc + nvirt - 1);
  Mfull.Fill(aux, dft, orb.MOs().eigenvectors());
  RPA rpa_full(log, Mfull);
  rpa_full.configure(nocc - 1, 0, nocc + nvirt - 1);
  rpa_full.setRPAInputEnergies(orb.MOs().eigenvalues());

  // embedded levels without the inactive ones shifted to the top
  const Eigen::MatrixXd Cemb = emb.eigenvectors().leftCols(nA + nvirt);
  TCMatrix_gwbse Memb;
  Memb.Initialize(aux.AOBasisSize(), 0, nA + nvirt - 1, 0, nA + nvirt - 1);
  Memb.Fill(aux, dft, Cemb);
  RPA rpa_emb(log, Memb);
  rpa_emb.configure(nA - 1, 0, nA + nvirt - 1);
  rpa_emb.setRPAInputEnergies(emb.eigenvalues().head(nA + nvirt));

  const InactiveScreeningInput in =
      InactiveScreening::Prepare(orb, InactiveScreeningInput::Mode::Exact);
  const double eta = rpa_emb.getEta();
  const double wi = 0.5;
  Memb.setAdditionalResponse(
      {{false, 0.0,
        InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9, eta)},
       {true, wi,
        InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9, wi)}});
  const Eigen::MatrixXd e0 = rpa_full.calculate_epsilon_r(0.0);
  const Eigen::MatrixXd ei = rpa_full.calculate_epsilon_i(wi);
  BOOST_CHECK_SMALL(MaxDiff(rpa_emb.calculate_epsilon_r(0.0), e0), 1e-10);
  BOOST_CHECK_SMALL(MaxDiff(rpa_emb.calculate_epsilon_i(wi), ei), 1e-10);

  // rotated auxiliary index
  const Index naux = e0.rows();
  const Eigen::MatrixXd U =
      Eigen::HouseholderQR<Eigen::MatrixXd>(Eigen::MatrixXd::Random(naux, naux))
          .householderQ();
  Memb.MultiplyRightWithAuxMatrix(U);
  BOOST_CHECK_SMALL(
      MaxDiff(rpa_emb.calculate_epsilon_i(wi), U.transpose() * ei * U), 1e-10);

  // other frequencies have no inactive response
  BOOST_CHECK_THROW(rpa_emb.calculate_epsilon_i(0.3), std::runtime_error);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(inactive_poles_in_the_rpa_pair_problem) {
  // block mode: the inactive response fitted by poles at 0 and 0.5i, which
  // enter the RPA two-particle Hamiltonian (Sigma_Exact) as
  // pseudo-transitions
  libint2::initialize();
  Orbitals orb = EmbeddedWater();
  const AOBasis& dft = orb.getDftBasis();
  const AOBasis aux = AuxBasis(orb);
  const Index nocc = orb.getNumberOfAlphaElectrons();
  const Index nvirt = orb.MOs().eigenvectors().cols() - nocc;
  const Index nA = orb.getNumOfActiveElectrons() / 2;
  const tools::EigenSystem& emb = orb.getEmbeddedMOs();
  Logger log;

  const InactiveScreeningInput in =
      InactiveScreening::Prepare(orb, InactiveScreeningInput::Mode::Block);
  const double wi = 0.5;
  const Eigen::MatrixXd Y0 =
      InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9, 0.0);
  const Eigen::MatrixXd Yi =
      InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9, wi);
  const InactiveScreening::PoleFit fit =
      InactiveScreening::FitPoles(Y0, Yi, wi);
  BOOST_CHECK_EQUAL(fit.clipped, 0);
  BOOST_TEST_MESSAGE("block: " << fit.energies.size() << " poles, "
                               << fit.dropped << " directions dropped");
  // the poles reproduce both points (on the range of Y(0))
  auto poles_at = [&](double w) {
    const Eigen::ArrayXd W = fit.energies.array();
    const Eigen::VectorXd g = 4 * W / (W.square() + w * w);
    return Eigen::MatrixXd(fit.vectors * g.asDiagonal() *
                           fit.vectors.transpose());
  };
  BOOST_CHECK_SMALL(MaxDiff(poles_at(0.0), Y0) / Y0.cwiseAbs().maxCoeff(),
                    1e-8);
  BOOST_CHECK_SMALL(MaxDiff(poles_at(wi), Yi) / Yi.cwiseAbs().maxCoeff(), 1e-8);
  // between the points the two-point fit is an approximation
  const Eigen::MatrixXd Y3 =
      InactiveScreening::Kernel(aux, dft, in, emb, nA - 1, 4e9, 0.3);
  BOOST_TEST_MESSAGE("relative |Y_poles(0.3i) - Y(0.3i)| = "
                     << MaxDiff(poles_at(0.3), Y3) / Y3.cwiseAbs().maxCoeff());

  // embedded active integrals with the poles
  const Eigen::MatrixXd Cemb = emb.eigenvectors().leftCols(nA + nvirt);
  TCMatrix_gwbse Memb;
  Memb.Initialize(aux.AOBasisSize(), 0, nA + nvirt - 1, 0, nA + nvirt - 1);
  Memb.Fill(aux, dft, Cemb);
  Memb.setAdditionalPoles(fit.energies, fit.vectors);
  RPA rpa(log, Memb);
  rpa.configure(nA - 1, 0, nA + nvirt - 1);
  rpa.setRPAInputEnergies(emb.eigenvalues().head(nA + nvirt));
  const Index naux = Y0.rows();
  const Eigen::MatrixXd one = Eigen::MatrixXd::Identity(naux, naux);

  // epsilon with the poles at any frequency
  Memb.setAdditionalPoles(Eigen::VectorXd(), Eigen::MatrixXd());
  const Eigen::MatrixXd eps_A = rpa.calculate_epsilon_i(0.3);
  Memb.setAdditionalPoles(fit.energies, fit.vectors);
  BOOST_CHECK_SMALL(
      MaxDiff(rpa.calculate_epsilon_i(0.3), eps_A + poles_at(0.3)), 1e-12);

  // the pair problem with pseudo-transitions gives the poles of
  // eps^-1 - 1 = -4 Z diag(W/(nu^2 + W^2)) Z^T, also in a rotated frame
  auto check_pair_problem = [&]() {
    const RPA::rpa_eigensolution sol = rpa.Diagonalize_H2p();
    const Index n_unocc = nvirt;
    const Index rpasize = nA * n_unocc;
    BOOST_REQUIRE_EQUAL(sol.omega.size(), rpasize + fit.energies.size());
    Eigen::MatrixXd M(sol.omega.size(), naux);
    for (Index v = 0; v < nA; v++) {
      M.middleRows(v * n_unocc, n_unocc) = Memb[v].middleRows(nA, n_unocc);
    }
    M.bottomRows(fit.energies.size()) =
        Memb.AdditionalPoleVectors().transpose();
    const Eigen::MatrixXd Z = M.transpose() * sol.XpY;
    for (double nu : {0.0, 0.3, 1.7}) {
      const Eigen::ArrayXd O = sol.omega.array();
      const Eigen::VectorXd g = -4 * O / (O.square() + nu * nu);
      const Eigen::MatrixXd from_poles = Z * g.asDiagonal() * Z.transpose();
      const Eigen::MatrixXd direct =
          InverseSPD(rpa.calculate_epsilon_i(nu)) - one;
      BOOST_CHECK_SMALL(MaxDiff(from_poles, direct), 1e-8);
    }
  };
  check_pair_problem();
  const Eigen::MatrixXd U =
      Eigen::HouseholderQR<Eigen::MatrixXd>(Eigen::MatrixXd::Random(naux, naux))
          .householderQ();
  Memb.MultiplyRightWithAuxMatrix(U);
  check_pair_problem();
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(combined_with_environment) {
  // 1 + R = ((1 + R_env)^-1 + X_B)^-1 for a synthetic environment
  const Index n = 6;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd X_B = 0.3 * A * A.transpose();
  Eigen::MatrixXd Bm = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd R_env = -0.05 * Bm * Bm.transpose();
  const Eigen::MatrixXd one = Eigen::MatrixXd::Identity(n, n);
  const Eigen::MatrixXd R = InactiveScreening::ReactionField(X_B, R_env);
  const Eigen::MatrixXd lhs = (one + R).inverse();
  const Eigen::MatrixXd rhs = (one + R_env).inverse() + X_B;
  BOOST_CHECK_SMALL(MaxDiff(lhs, rhs), 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
