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

#define BOOST_TEST_MODULE ppm_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/ppm.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(ppm_test)

BOOST_AUTO_TEST_CASE(ppm_full) {
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/ppm/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/ppm/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Orbitals orb;
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);

  AOKinetic kinetic;
  kinetic.Fill(aobasis);

  AOMultipole esp;
  esp.FillPotential(aobasis, orbitals.QMAtoms());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kinetic.Matrix() +
                                                    esp.Matrix());

  Logger log;
  TCMatrix_gwbse Mmn{log};
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, es.eigenvectors());

  RPA rpa = RPA(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(es.eigenvalues());

  PPM ppm;
  ppm.PPM_construct_parameters(rpa);

  Eigen::VectorXd ppm_freq = Eigen::VectorXd::Zero(17);
  ppm_freq << 19.4503, 12.2429, 10.2167, 10.2167, 10.2167, 17.0832, 12.7117,
      12.7117, 12.7117, 10.9455, 10.9455, 10.9455, 11.1861, 9.60843, 9.60843,
      9.60843, 9.61518;
  Eigen::VectorXd ppm_w = Eigen::VectorXd::Zero(17);
  ppm_w << 3.59422e-05, 0.00121795, 0.00343632, 0.00343632, 0.00343632,
      0.00394659, 0.012066, 0.012066, 0.012066, 0.0241118, 0.0241118, 0.0241118,
      0.0286786, 0.11036, 0.11036, 0.11036, 0.191732;

  bool f_check = ppm_freq.isApprox(ppm.getPpm_freq(), 0.0001);

  bool w_check = ppm_w.isApprox(ppm.getPpm_weight(), 0.0001);

  if (!f_check) {
    cout << "ppm_freq" << endl;
    cout << ppm.getPpm_freq() << endl;
    cout << "ppm_freq_ref" << endl;
    cout << ppm_freq << endl;
  }

  if (!w_check) {
    cout << "ppm_w" << endl;
    cout << ppm.getPpm_weight() << endl;
    cout << "ppm_w_ref" << endl;
    cout << ppm_w << endl;
  }

  BOOST_CHECK_EQUAL(f_check, 1);
  BOOST_CHECK_EQUAL(w_check, 1);
}

BOOST_AUTO_TEST_SUITE_END()
