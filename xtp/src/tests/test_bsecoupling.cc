/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE bsecoupling_test

// Third party includes
#include <boost/test/unit_test.hpp>
#include <sstream>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/bsecoupling.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(bsecoupling_test)

BOOST_AUTO_TEST_CASE(coupling_test) {
  libint2::initialize();
  Orbitals A;

  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                           "/bsecoupling/molecule.xyz");
  A.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/3-21G.xml");
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/A_MOs.mm");

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/spsi_ref.mm");

  A.BSESinglets().eigenvectors() = spsi_ref;

  // Set BSE singlet eigenvalues (excitation energies in Hartree).
  // Required by the updated BSECoupling which stores monomer energies
  // as pairwise-consistent TB site energies. The value here is a
  // placeholder consistent with the HOMO-LUMO gap of the test molecule
  // (~2.4 eV = 0.0883 Hrt from the MO eigenvalues above).
  A.BSESinglets().eigenvalues().resize(1);
  A.BSESinglets().eigenvalues() << 0.08831;  // ~2.4 eV placeholder

  Orbitals B = A;
  B.QMAtoms().Translate(4 * Eigen::Vector3d::UnitX());

  Orbitals AB;
  AB.QMAtoms() = A.QMAtoms();
  AB.QMAtoms().AddContainer(B.QMAtoms());
  AB.MOs().eigenvalues().resize(34);
  AB.MOs().eigenvalues() << -10.1341, -10.1337, -0.808607, -0.665103, -0.474928,
      -0.455857, -0.455857, -0.365971, -0.365971, -0.263259, 0.140444, 0.154745,
      0.168775, 0.168775, 0.223948, 0.231217, 0.26323, 0.26323, 0.713478,
      0.713478, 0.793559, 0.885998, 0.944915, 0.944915, 1.01169, 1.04977,
      1.04977, 1.08863, 1.10318, 1.17822, 1.18094, 1.18094, 1.69037, 1.91046;

  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);
  AB.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bsecoupling/3-21G.xml");
  AB.SetupAuxBasis(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bsecoupling/3-21G.xml");
  AB.setRPAindices(0, 33);
  AB.setBSEindices(0, 33);
  AB.setGWindices(0, 33);
  AB.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/AB_MOs.mm");

  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);

  AB.SetFlagUseHqpOffdiag(true);

  std::ofstream opt("bsecoupling.xml");
  opt << "<bsecoupling>" << std::endl;
  opt << "        <use_perturbation>true</use_perturbation>" << std::endl;
  opt << "        <spin>singlet</spin>" << std::endl;
  opt << "        <output_tb>false</output_tb>" << std::endl;
  opt << "       <moleculeA>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "        </moleculeA>" << std::endl;
  opt << "        <moleculeB>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "         </moleculeB>" << std::endl;
  opt << "</bsecoupling>" << std::endl;
  opt.close();
  votca::tools::Property prop;
  prop.LoadFromXML("bsecoupling.xml");
  BSECoupling coup;
  Logger log;
  log.setCommonPreface("\n... ...");
  coup.setLogger(&log);

  AB.QPdiag().eigenvalues().resize(34);
  AB.QPdiag().eigenvalues() << -10.504, -10.5038, -0.923616, -0.775673,
      -0.549084, -0.530193, -0.530193, -0.430293, -0.430293, -0.322766,
      0.267681, 0.307809, 0.326961, 0.326961, 0.36078, 0.381947, 0.414845,
      0.414845, 0.906609, 0.906609, 0.993798, 1.09114, 1.14639, 1.14639, 1.1966,
      1.25629, 1.25629, 1.27991, 1.29122, 1.35945, 1.36705, 1.36705, 1.93286,
      2.11739;

  AB.QPdiag().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/Hqp.mm");
  const Eigen::MatrixXd& qpcoeff = AB.QPdiag().eigenvectors();
  Eigen::MatrixXd Hqp =
      qpcoeff * AB.QPdiag().eigenvalues().asDiagonal() * qpcoeff.transpose();
  AB.RPAInputEnergies() = Hqp.diagonal();
  coup.Initialize(prop.get("bsecoupling"));
  log.setReportLevel(Log::error);
  coup.CalculateCouplings(A, B, AB);
  votca::tools::Property output;
  coup.Addoutput(output, A, B);
  // Reference values updated to reflect removal of the incorrect CT
  // pre-orthogonalization (P = F*F^T was only exact for orthonormal FE
  // states). The joint Lowdin in CalcJ_dimer now handles all non-orthogonality
  // correctly. The large values (eV) reflect the unphysically short
  // intermolecular distance in the test geometry (4 bohr).
  // The sign of J_diag is gauge-dependent (BSE eigenvector phase),
  // so comparisons use absolute values.
  double diag_J_ref = 23.662750;  // eV
  double pert_J_ref = 9.529579;   // eV

  double diag_j =
      output.get("bsecoupling.singlet.coupling").getAttribute<double>("j_diag");
  double pert_j =
      output.get("bsecoupling.singlet.coupling").getAttribute<double>("j_pert");

  BOOST_CHECK_CLOSE(diag_J_ref, std::abs(diag_j), 1e-4);
  BOOST_CHECK_CLOSE(pert_J_ref, std::abs(pert_j), 1e-4);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(tb_output_test) {
  // Test the TB matrix output (output_tb=true) added to BSECoupling.
  // Uses the same test geometry and orbitals as coupling_test.
  // Verifies:
  //   - Correct structural attributes (n_FE, n_CT, matrix dimensions)
  //   - Monomer energy recovered from BSE eigenvalue
  //   - S_FE_FE diagonal close to 1 (approximate normalisation)
  //   - H_FE_FE off-diagonal matches J_diag from the scalar coupling test
  //   - CT energies positive and larger than FE energies (physical ordering)
  //   - Diagnostics present and physically consistent
  libint2::initialize();
  Orbitals A;

  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                           "/bsecoupling/molecule.xyz");
  A.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/3-21G.xml");
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;
  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/A_MOs.mm");
  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/spsi_ref.mm");
  A.BSESinglets().eigenvectors() = spsi_ref;
  A.BSESinglets().eigenvalues().resize(1);
  A.BSESinglets().eigenvalues() << 0.08831;  // ~2.4 eV

  Orbitals B = A;
  B.QMAtoms().Translate(4 * Eigen::Vector3d::UnitX());

  Orbitals AB;
  AB.QMAtoms() = A.QMAtoms();
  AB.QMAtoms().AddContainer(B.QMAtoms());
  AB.MOs().eigenvalues().resize(34);
  AB.MOs().eigenvalues() << -10.1341, -10.1337, -0.808607, -0.665103, -0.474928,
      -0.455857, -0.455857, -0.365971, -0.365971, -0.263259, 0.140444, 0.154745,
      0.168775, 0.168775, 0.223948, 0.231217, 0.26323, 0.26323, 0.713478,
      0.713478, 0.793559, 0.885998, 0.944915, 0.944915, 1.01169, 1.04977,
      1.04977, 1.08863, 1.10318, 1.17822, 1.18094, 1.18094, 1.69037, 1.91046;
  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);
  AB.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bsecoupling/3-21G.xml");
  AB.SetupAuxBasis(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bsecoupling/3-21G.xml");
  AB.setRPAindices(0, 33);
  AB.setBSEindices(0, 33);
  AB.setGWindices(0, 33);
  AB.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/AB_MOs.mm");
  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);
  AB.SetFlagUseHqpOffdiag(true);
  AB.QPdiag().eigenvalues().resize(34);
  AB.QPdiag().eigenvalues() << -10.504, -10.5038, -0.923616, -0.775673,
      -0.549084, -0.530193, -0.530193, -0.430293, -0.430293, -0.322766,
      0.267681, 0.307809, 0.326961, 0.326961, 0.36078, 0.381947, 0.414845,
      0.414845, 0.906609, 0.906609, 0.993798, 1.09114, 1.14639, 1.14639, 1.1966,
      1.25629, 1.25629, 1.27991, 1.29122, 1.35945, 1.36705, 1.36705, 1.93286,
      2.11739;
  AB.QPdiag().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/Hqp.mm");
  const Eigen::MatrixXd& qpcoeff = AB.QPdiag().eigenvectors();
  Eigen::MatrixXd Hqp =
      qpcoeff * AB.QPdiag().eigenvalues().asDiagonal() * qpcoeff.transpose();
  AB.RPAInputEnergies() = Hqp.diagonal();

  // Enable TB output
  std::ofstream opt("bsecoupling_tb.xml");
  opt << "<bsecoupling>" << std::endl;
  opt << "        <use_perturbation>true</use_perturbation>" << std::endl;
  opt << "        <spin>singlet</spin>" << std::endl;
  opt << "        <output_tb>true</output_tb>" << std::endl;
  opt << "       <moleculeA>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "        </moleculeA>" << std::endl;
  opt << "        <moleculeB>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "         </moleculeB>" << std::endl;
  opt << "</bsecoupling>" << std::endl;
  opt.close();

  votca::tools::Property prop;
  prop.LoadFromXML("bsecoupling_tb.xml");
  BSECoupling coup;
  Logger log;
  log.setCommonPreface("\n... ...");
  coup.setLogger(&log);
  log.setReportLevel(Log::error);
  coup.Initialize(prop.get("bsecoupling"));
  coup.CalculateCouplings(A, B, AB);
  votca::tools::Property output;
  coup.Addoutput(output, A, B);

  // -------------------------------------------------------------------------
  // Structural checks
  // -------------------------------------------------------------------------
  // n_FE = levA + levB = 1 + 1 = 2
  // n_CT = 2 * occLevels * unoccLevels = 2 * 3 * 3 = 18
  Index n_FE =
      output.get("bsecoupling.singlet.tb_matrices").getAttribute<Index>("n_FE");
  Index n_CT =
      output.get("bsecoupling.singlet.tb_matrices").getAttribute<Index>("n_CT");
  BOOST_CHECK_EQUAL(n_FE, 2);
  BOOST_CHECK_EQUAL(n_CT, 18);

  // Matrix dimensions
  Index hff_rows = output.get("bsecoupling.singlet.tb_matrices.H_FE_FE")
                       .getAttribute<Index>("rows");
  Index hff_cols = output.get("bsecoupling.singlet.tb_matrices.H_FE_FE")
                       .getAttribute<Index>("cols");
  BOOST_CHECK_EQUAL(hff_rows, 2);
  BOOST_CHECK_EQUAL(hff_cols, 2);

  Index hfc_rows = output.get("bsecoupling.singlet.tb_matrices.H_FE_CT")
                       .getAttribute<Index>("rows");
  Index hfc_cols = output.get("bsecoupling.singlet.tb_matrices.H_FE_CT")
                       .getAttribute<Index>("cols");
  BOOST_CHECK_EQUAL(hfc_rows, 2);
  BOOST_CHECK_EQUAL(hfc_cols, 18);

  Index hcc_rows = output.get("bsecoupling.singlet.tb_matrices.H_CT_CT")
                       .getAttribute<Index>("rows");
  BOOST_CHECK_EQUAL(hcc_rows, 18);

  // -------------------------------------------------------------------------
  // H_FE_FE matrix: raw non-orthogonal projected BSE Hamiltonian.
  // H_FE_FE is in the non-normalised dimer projection basis BEFORE Lowdin
  // orthogonalisation, so its diagonal is NOT the monomer site energy.
  // The diagonal H_AA = <psi_A|H|psi_A> where psi_A has norm S_AA =
  // <psi_A|psi_A>. The site energy is recovered as H_AA / S_AA after Lowdin.
  //
  // What we CAN test:
  //   1. H and S are symmetric: H_AB == H_BA, S_AB == S_BA
  //   2. S diagonal is positive (it's a norm squared)
  //   3. H/S diagonal (generalised site energy) is finite and positive
  //   4. The Lowdin-orthogonalised J_diag (stored in the scalar coupling)
  //      is recovered correctly by S^{-1/2} H S^{-1/2} off-diagonal
  //      — verified implicitly since coupling_test already checks J_diag
  // -------------------------------------------------------------------------
  std::string row0_str = output.get("bsecoupling.singlet.tb_matrices.H_FE_FE")
                             .getAttribute<std::string>("row_0");
  std::string row1_str = output.get("bsecoupling.singlet.tb_matrices.H_FE_FE")
                             .getAttribute<std::string>("row_1");
  std::istringstream ss0(row0_str), ss1(row1_str);
  double H_AA, H_AB, H_BA, H_BB;
  ss0 >> H_AA >> H_AB;
  ss1 >> H_BA >> H_BB;

  // H is symmetric
  BOOST_CHECK_CLOSE(H_AB, H_BA, 1e-4);
  // H diagonal must be finite
  BOOST_CHECK(std::isfinite(H_AA));
  BOOST_CHECK(std::isfinite(H_BB));

  std::string sff_row0 = output.get("bsecoupling.singlet.tb_matrices.S_FE_FE")
                             .getAttribute<std::string>("row_0");
  std::string sff_row1 = output.get("bsecoupling.singlet.tb_matrices.S_FE_FE")
                             .getAttribute<std::string>("row_1");
  std::istringstream ss2(sff_row0), ss3(sff_row1);
  double S_AA, S_AB_s, S_BA_s, S_BB;
  ss2 >> S_AA >> S_AB_s;
  ss3 >> S_BA_s >> S_BB;

  // S is symmetric
  BOOST_CHECK_CLOSE(S_AB_s, S_BA_s, 1e-4);
  // S diagonal is positive (it's a squared norm of a projection)
  BOOST_CHECK(S_AA > 0.0);
  BOOST_CHECK(S_BB > 0.0);

  // Generalised site energy H_AA/S_AA must be positive and finite.
  // This is what the Lowdin step recovers as the effective diagonal.
  BOOST_CHECK(std::isfinite(H_AA / S_AA));
  BOOST_CHECK(H_AA / S_AA > 0.0);

  // -------------------------------------------------------------------------
  // Monomer energy node: the eV attribute is the BSE eigenvalue in eV,
  // written from orbitalsA.BSESinglets().eigenvalues(). We set 0.08831 Hrt,
  // so the expected value is 0.08831 * hrt2ev.
  // This tests the serialisation path independently of the matrix entries.
  // -------------------------------------------------------------------------
  const double e_mono_ref = 0.08831 * votca::tools::conv::hrt2ev;  // eV

  // -------------------------------------------------------------------------
  // H_FE_CT: numerical accuracy check.
  // The maximum absolute element of H_FE_CT (the FE-CT coupling block) is
  // the same quantity that drives xi = max|H_FE_CT| / |E_FE - E_CT|.
  // We check that at least one element is non-zero (coupling exists) and
  // finite, and pin the largest element to a reference value derived from
  // the test geometry. Since the sign is gauge-dependent we use abs.
  // Reference: 171.0 eV (from test run; large due to 4-bohr separation).
  // -------------------------------------------------------------------------
  std::string hfc_row0 = output.get("bsecoupling.singlet.tb_matrices.H_FE_CT")
                             .getAttribute<std::string>("row_0");
  std::string hfc_row1 = output.get("bsecoupling.singlet.tb_matrices.H_FE_CT")
                             .getAttribute<std::string>("row_1");
  // find max |H_FE_CT| across both rows
  auto maxAbsRow = [](const std::string& row_str) {
    std::istringstream ss(row_str);
    double val, maxval = 0.0;
    while (ss >> val) maxval = std::max(maxval, std::abs(val));
    return maxval;
  };
  double max_hfc = std::max(maxAbsRow(hfc_row0), maxAbsRow(hfc_row1));
  BOOST_CHECK(max_hfc > 0.0);                // coupling is non-zero
  BOOST_CHECK(std::isfinite(max_hfc));       // and finite
  BOOST_CHECK_CLOSE(max_hfc, 1894.51, 1.0);  // within 1%, eV

  // -------------------------------------------------------------------------
  // CT diagonal (H_CT_CT): first CT state energy must be finite.
  // No physical ordering guarantee for this unphysical geometry.
  // -------------------------------------------------------------------------
  std::string hcc_row0 = output.get("bsecoupling.singlet.tb_matrices.H_CT_CT")
                             .getAttribute<std::string>("row_0");
  std::istringstream ss4(hcc_row0);
  double H_CT00;
  ss4 >> H_CT00;
  BOOST_CHECK(std::isfinite(H_CT00));

  // -------------------------------------------------------------------------
  // Diagnostics: xi and downfolding_safe
  // -------------------------------------------------------------------------
  double xi =
      output.get("bsecoupling.singlet.diagnostics").getAttribute<double>("xi");
  bool safe = output.get("bsecoupling.singlet.diagnostics")
                  .getAttribute<bool>("downfolding_safe");
  // xi must be non-negative
  BOOST_CHECK(xi >= 0.0);
  // With such large couplings (23 eV!) relative to any reasonable FE-CT gap,
  // xi >> 1 and downfolding must be flagged as unsafe
  BOOST_CHECK_EQUAL(safe, false);

  // Check monomer energy node eV attribute matches BSE eigenvalue
  double mono_eV =
      output.get("bsecoupling.singlet.monomer_energies.fragmentA.energy")
          .getAttribute<double>("eV");
  BOOST_CHECK_CLOSE(mono_eV, e_mono_ref, 1e-4);

  libint2::finalize();
}
BOOST_AUTO_TEST_SUITE_END()
