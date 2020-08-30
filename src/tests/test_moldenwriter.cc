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

#define BOOST_TEST_MODULE moldenwriter_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/moldenreader.h"
#include "votca/xtp/moldenwriter.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(moldenwriter_test)

BOOST_AUTO_TEST_CASE(moldenwriter_test) {

  // Setup orbitals object
  Orbitals orbitals_ref;
  orbitals_ref.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                      "/molden/benzene.xyz");
  orbitals_ref.MOs().eigenvectors() =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/molden/orbitalsMOs_ref.mm");
  orbitals_ref.setBasisSetSize(222);
  orbitals_ref.setNumberOfOccupiedLevels(21);
  Eigen::VectorXd eigenvalues(222);
  eigenvalues << -10.2307, -10.2305, -10.2304, -10.2299, -10.2298, -10.2296,
      -0.879257, -0.768798, -0.768711, -0.621806, -0.621652, -0.540065,
      -0.476447, -0.46157, -0.434835, -0.434726, -0.383707, -0.356208,
      -0.356125, -0.267169, -0.267002, -0.00952, -0.008996, 0.05598, 0.088658,
      0.090001, 0.114775, 0.11618, 0.138134, 0.147314, 0.160622, 0.218345,
      0.222561, 0.223444, 0.226297, 0.242926, 0.243833, 0.259218, 0.259671,
      0.274361, 0.274586, 0.331821, 0.337298, 0.337413, 0.363761, 0.36524,
      0.388862, 0.410777, 0.423401, 0.423908, 0.442365, 0.442973, 0.447599,
      0.464177, 0.474317, 0.494663, 0.498432, 0.499962, 0.515936, 0.516765,
      0.65182, 0.653069, 0.654124, 0.689895, 0.692011, 0.74763, 0.748273,
      0.767789, 0.768662, 0.800673, 0.844666, 0.844947, 0.858436, 0.895223,
      0.904508, 1.01265, 1.01384, 1.0207, 1.02341, 1.05627, 1.05669, 1.17412,
      1.18522, 1.18744, 1.21225, 1.21405, 1.22555, 1.22704, 1.25405, 1.25509,
      1.29422, 1.386, 1.38725, 1.41053, 1.46334, 1.46471, 1.49732, 1.49796,
      1.5017, 1.50212, 1.60217, 1.61954, 1.6808, 1.68477, 1.68609, 1.68685,
      1.69742, 1.69751, 1.71942, 1.72043, 1.77277, 1.85733, 1.90792, 1.90824,
      2.03662, 2.03678, 2.07245, 2.08436, 2.09265, 2.09278, 2.16536, 2.16553,
      2.23367, 2.24925, 2.24947, 2.32663, 2.45893, 2.50001, 2.50218, 2.56338,
      2.56399, 2.56644, 2.57395, 2.59273, 2.65586, 2.65634, 2.67822, 2.67849,
      2.67853, 2.68061, 2.68102, 2.7492, 2.75594, 2.75676, 2.80601, 2.8122,
      2.81505, 2.85033, 2.85264, 2.85362, 2.85602, 2.85902, 2.9094, 3.02012,
      3.10981, 3.11107, 3.11211, 3.11269, 3.139, 3.14012, 3.14726, 3.14827,
      3.20044, 3.22395, 3.22432, 3.25389, 3.25987, 3.26037, 3.32212, 3.36218,
      3.37524, 3.3753, 3.39566, 3.45773, 3.45849, 3.46522, 3.50579, 3.52481,
      3.52488, 3.5959, 3.65255, 3.65326, 3.75155, 3.7529, 3.75758, 3.80608,
      3.88474, 4.00529, 4.00563, 4.05716, 4.07444, 4.0758, 4.1153, 4.19871,
      4.19913, 4.21241, 4.21886, 4.21939, 4.27827, 4.27872, 4.29987, 4.30042,
      4.62163, 4.64957, 4.65574, 4.65786, 4.83063, 4.83152, 4.89966, 4.90036,
      5.19347, 5.19352, 5.28327, 5.30115, 5.41754, 5.63135, 21.9515, 22.5924,
      22.5926, 22.6859, 22.6865, 23.1019;
  orbitals_ref.MOs().eigenvalues() = eigenvalues;
  orbitals_ref.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                               "/molden/def2-tzvp.xml");
  orbitals_ref.setAuxbasisName("aux-def2-tzvp");

  // write orbitals object to molden file
  Logger log;
  MoldenWriter moldenWriter(log);
  moldenWriter.WriteFile("moldenFile.molden", orbitals_ref);

  // read in written molden file
  MoldenReader molden(log);
  molden.setBasissetInfo(
      std::string(XTP_TEST_DATA_FOLDER) + "/molden/def2-tzvp.xml",
      "aux-def2-tzvp");
  Orbitals orbitals;
  molden.parseMoldenFile("moldenFile.molden", orbitals);

  // Check if MO's are equal
  BOOST_CHECK(orbitals_ref.MOs().eigenvectors().isApprox(
      orbitals.MOs().eigenvectors(), 1e-5));

  // Check if atoms are equal
  BOOST_CHECK(orbitals.QMAtoms().size() == orbitals_ref.QMAtoms().size());
  for (int i = 0; i < orbitals.QMAtoms().size(); i++) {
    BOOST_CHECK(orbitals.QMAtoms()[i].getPos().isApprox(
        orbitals_ref.QMAtoms()[i].getPos(), 1e-3));
  }
}
}