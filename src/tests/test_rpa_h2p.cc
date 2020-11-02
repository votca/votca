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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE rpa_test

// Third party includes
#include <boost/test/unit_test.hpp>
// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"
#include <libint2/initialize.h>
using namespace std;
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_h2p) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/rpa/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/rpa/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd eigenvals = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvals.mm");

  Eigen::MatrixXd eigenvectors = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvectors.mm");

  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, eigenvectors);

  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(eigenvals);
  rpa.configure(4, 0, 16);

  Eigen::VectorXd rpa_omega_ref = Eigen::VectorXd(60);
  rpa_omega_ref << 0.104192, 0.104192, 0.187814, 0.559693, 0.559693, 0.572575,
      0.577988, 0.577989, 0.579088, 0.618403, 0.618403, 0.67005, 0.678538,
      0.678538, 0.722771, 1.10797, 1.41413, 1.41413, 1.58866, 1.60381, 1.60381,
      1.64709, 1.87331, 1.87331, 1.88646, 1.8926, 1.89268, 1.89268, 1.933,
      1.933, 2.01832, 2.40974, 2.42192, 2.42192, 2.46371, 2.85829, 2.8853,
      2.8853, 2.90367, 2.90367, 2.92541, 2.94702, 3.3382, 3.3382, 3.35102,
      3.3566, 3.3566, 3.35835, 3.39617, 3.39617, 4.22882, 4.71607, 4.72233,
      4.72233, 4.76567, 16.5917, 17.0793, 17.093, 17.093, 17.1377;

  Eigen::VectorXd rpa_XpY_diag_ref = Eigen::VectorXd(60);
  rpa_XpY_diag_ref << 0.00898067, -0.00898068, 0.00228621, -7.40059e-10,
      -0.000447591, -7.21786e-10, 2.86571e-09, 1.16917e-08, 5.18349e-09,
      0.000124051, 1.89797e-10, -1.549e-05, -0.00513947, -0.00513956,
      -7.71176e-08, -2.09605e-08, 0.00412762, 0.00412738, 6.253e-09,
      3.93273e-05, 2.09406e-05, -0.000532391, 3.70347e-05, 2.0857e-06,
      -6.13653e-10, 0.000723575, 3.41096e-05, 0.00942648, 0.0497694, 0.0497693,
      -1.28726e-07, 3.68417e-09, -7.00924e-06, 7.01603e-06, -1.36252e-09,
      -1.08353e-09, 0.000549202, 0.000549203, -6.22315e-09, -3.83321e-08,
      -2.4032e-08, -2.19679e-09, 5.76794e-10, -1.46987e-08, -1.27183e-08,
      -0.138355, -6.49432e-09, -2.71175e-05, 2.69024e-05, -2.69024e-05,
      0.000155952, 0.000311876, -0.0002498, -0.000249804, 0.000816289,
      1.05559e-06, 1.84996e-12, 5.76753e-06, -9.90455e-11, -0.000800706;

  RPA::rpa_eigensolution sol = rpa.Diagonalize_H2p();

  BOOST_CHECK_CLOSE(sol.ERPA_correlation, -0.0587973, 1e-4);

  bool check_rpa_eigenvalues = rpa_omega_ref.isApprox(sol.omega, 0.0001);
  if (!check_rpa_eigenvalues) {
    cout << "rpa_omega" << endl;
    cout << sol.omega << endl;
    cout << "rpa_omega_ref" << endl;
    cout << rpa_omega_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_rpa_eigenvalues, 1);

  Eigen::VectorXd rpa_XpY_diag = sol.XpY.diagonal();

  bool check_rpa_XpY_diag =
      rpa_XpY_diag_ref.cwiseAbs().isApprox(rpa_XpY_diag.cwiseAbs(), 0.0001);
  if (!check_rpa_XpY_diag) {
    cout << "rpa_XpY_diag" << endl;
    cout << rpa_XpY_diag.transpose() << endl;
    cout << "rpa_XpY_diag_ref" << endl;
    cout << rpa_XpY_diag_ref.transpose() << endl;
  }
  BOOST_CHECK_EQUAL(check_rpa_XpY_diag, 1);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
