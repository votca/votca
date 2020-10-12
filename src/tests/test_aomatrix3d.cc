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

#define BOOST_TEST_MODULE aomatrix3d_test

// Third party includes
#include <boost/test/unit_test.hpp>

#include <votca/tools/eigenio_matrixmarket.h>
// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/orbitals.h"
#include <libint2/initialize.h>
using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aomatrix3d_test)

BOOST_AUTO_TEST_CASE(aomatrices3d_test) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aomatrix3d/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  AODipole dip;
  dip.Fill(aobasis);
  std::array<Eigen::MatrixXd, 3> dip_ref;
  for (unsigned i = 0; i < dip_ref.size(); i++) {
    dip_ref[i] = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/dip_ref" +
        std::to_string(i) + ".mm");
  }

  for (unsigned i = 0; i < dip_ref.size(); i++) {
    bool check_dip = dip.Matrix()[i].isApprox(dip_ref[i], 0.0001);
    BOOST_CHECK_EQUAL(check_dip, true);
    if (!check_dip) {
      cout << "ref" << i << endl;
      cout << dip_ref[i] << endl;
      cout << "result" << i << endl;
      cout << dip.Matrix()[i] << endl;
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(large_l_test) {
  libint2::initialize();
  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/C2.xyz");

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/G.xml");
  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);

  AODipole dip;
  dip.Fill(dftbasis);

  std::array<Eigen::MatrixXd, 3> dip_ref;
  for (unsigned i = 0; i < dip_ref.size(); i++) {
    dip_ref[i] = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/dip_ref_large_" +
        std::to_string(i) + ".mm");
  }

  for (unsigned i = 0; i < dip_ref.size(); i++) {

    bool check_dip = dip.Matrix()[i].isApprox(dip_ref[i], 0.0001);
    BOOST_CHECK_EQUAL(check_dip, true);
    if (!check_dip) {
      cout << "ref" << i << endl;
      cout << dip_ref[i] << endl;
      cout << "result" << i << endl;
      cout << dip.Matrix()[i] << endl;
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
