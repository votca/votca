/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE density_filter_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/filterfactory.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(density_filter_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {
  FilterFactory::RegisterAll();
  std::unique_ptr<StateFilter_base> rho_f =
      std::unique_ptr<StateFilter_base>(Filter().Create("density"));

  std::ofstream opt("density_filter.xml");
  opt << "<density>0.0045</density>" << std::endl;
  opt.close();
  votca::tools::Property prop;
  prop.LoadFromXML("density_filter.xml");
  rho_f->Initialize(prop.get("density"));

  std::ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << std::endl;
  xyzfile << " methane" << std::endl;
  xyzfile << " C            .000000     .000000     .000000" << std::endl;
  xyzfile << " H            .629118     .629118     .629118" << std::endl;
  xyzfile << " H           -.629118    -.629118     .629118" << std::endl;
  xyzfile << " H            .629118    -.629118    -.629118" << std::endl;
  xyzfile << " H           -.629118     .629118    -.629118" << std::endl;
  xyzfile.close();

  std::ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << std::endl;
  basisfile << "  <element name=\"H\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "  <element name=\"C\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << std::endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell> " << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();

  Orbitals A;
  A.setDFTbasisName("3-21G.xml");
  A.QMAtoms().LoadFromFile("molecule.xyz");
  A.setBasisSetSize(17);
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.MOs().eigenvectors() = Eigen::MatrixXd::Zero(17, 17);
  A.MOs().eigenvectors() << -0.996559, -0.223082, 4.81443e-15, 2.21045e-15,
      -6.16146e-17, -3.16966e-16, 5.46703e-18, -1.09681e-15, -0.0301914,
      6.45993e-16, 1.05377e-16, 3.41154e-16, -0.102052, -5.79826e-16,
      9.38593e-16, -4.69346e-15, -0.111923, -0.0445146, 0.88316, -1.94941e-14,
      -8.67388e-15, -7.26679e-16, 1.16326e-14, -3.35886e-15, 2.37877e-14,
      0.866126, 3.2068e-15, 3.80914e-15, 3.24563e-15, -0.938329, -6.4404e-15,
      1.10811e-14, -5.5056e-14, -1.28767, 8.15798e-17, 2.30849e-14, 1.04169,
      0.117804, 0.0951759, 0.179467, 0.147031, 0.39183, -1.02927e-14, 0.32699,
      -0.223689, -0.130009, 1.0375e-15, -0.0940179, 0.126956, 0.0122904,
      1.41709e-15, 4.60157e-17, -7.1203e-15, 0.143338, -0.980459, -0.355251,
      0.41611, -0.10826, -0.149964, 2.41546e-16, 0.12214, -0.0512447, 0.39537,
      1.1054e-15, -0.0996828, -0.0636092, -0.105478, 5.10746e-15, -5.25872e-18,
      4.8424e-15, 0.0488925, 0.364515, -0.9863, 0.0447336, 0.417155, -0.177023,
      5.76117e-15, -0.228081, -0.348136, 0.0253377, -1.05286e-15, 0.079576,
      0.0703157, -0.117608, 5.31327e-15, 0.0472716, 0.235837, -3.58018e-15,
      -1.68354e-15, 2.3989e-15, -9.86879e-15, 4.52519e-15, -1.6106e-14,
      -0.599523, -1.31237e-14, -8.63443e-15, -8.61196e-15, 1.8333, 2.05275e-14,
      -3.9562e-14, 1.89874e-13, 4.24316, -2.74184e-16, -1.53939e-15, -0.162416,
      -0.0183675, -0.0148395, -0.151162, -0.123842, -0.330032, 1.10084e-15,
      -1.45092, 0.992556, 0.576875, -3.82954e-15, 0.604373, -0.816111,
      -0.0790061, -8.89474e-15, -2.24862e-16, 3.23655e-15, -0.0223487, 0.152869,
      0.0553894, -0.350483, 0.0911859, 0.126313, -5.48468e-15, -0.541961,
      0.227383, -1.75434, -3.89443e-15, 0.640788, 0.408897, 0.67804,
      -3.17156e-14, -2.81346e-17, -1.09423e-15, -0.00762313, -0.0568338,
      0.15378, -0.0376785, -0.351364, 0.149104, -4.94425e-15, 1.01204, 1.54475,
      -0.112429, 8.50653e-15, -0.511536, -0.452008, 0.756019, -3.3547e-14,
      -0.00106227, 0.0237672, 0.00345981, -0.00139675, -0.00349474, -0.597906,
      -0.425733, -0.0605479, -0.343823, 0.162103, -0.45692, 0.21318, -0.600309,
      0.310843, -0.36406, 0.574148, 0.0554949, -0.00928842, -0.0414346,
      0.0619999, -0.0250297, -0.0626259, 0.00227746, 0.00162164, 0.00023063,
      -0.0301047, 0.273177, -0.770004, 0.359253, 0.0095153, -0.8783, 1.02867,
      -1.62228, -1.24406, -0.00106227, 0.0237672, 0.00238182, 0.00205737,
      0.00402848, 0.262742, 0.151145, -0.671213, -0.343823, 0.317484, 0.12884,
      -0.40386, -0.600309, 0.201313, -0.327527, -0.641099, 0.0554949,
      -0.00928842, -0.0414346, 0.0426822, 0.0368682, 0.0721904, -0.0010008,
      -0.000575719, 0.00255669, -0.0301047, 0.535026, 0.217123, -0.680588,
      0.0095153, -0.568818, 0.925441, 1.81145, -1.24406, -0.00106227, 0.0237672,
      -0.00318563, 0.0034409, -0.00203628, 0.514364, -0.353326, 0.391148,
      -0.343823, -0.496623, -0.0536813, -0.176018, -0.600309, -0.744328,
      -0.01898, 0.0665156, 0.0554949, -0.00928842, -0.0414346, -0.0570866,
      0.0616609, -0.0364902, -0.00195924, 0.00134584, -0.0014899, -0.0301047,
      -0.836913, -0.0904642, -0.296627, 0.0095153, 2.10313, 0.0536287,
      -0.187943, -1.24406, -0.00106227, 0.0237672, -0.002656, -0.00410152,
      0.00150255, -0.1792, 0.627913, 0.340613, -0.343823, 0.0170366, 0.38176,
      0.366698, -0.600309, 0.232172, 0.710567, 0.000435528, 0.0554949,
      -0.00928842, -0.0414346, -0.0475955, -0.0734994, 0.0269257, 0.000682583,
      -0.00239176, -0.00129742, -0.0301047, 0.0287103, 0.643346, 0.617962,
      0.0095153, -0.656011, -2.00774, -0.0012306, -1.24406;

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = Eigen::MatrixXd::Zero(60, 3);
  spsi_ref << -0.00178316, -0.0558332, 0.0151767, 0.00568291, 0.0149417,
      0.0556358, 0.05758, -0.00320371, -0.00502105, 0.00379328, -0.00232255,
      -0.00817518, -0.00848959, -0.000618633, -0.00376334, -0.000395809,
      -0.00899117, 0.0023708, 7.13955e-08, 3.47498e-08, -1.08537e-08,
      0.000420413, 0.011896, -0.00320024, -0.00288487, 0.00320821, 0.0115465,
      0.0119767, 0.000355172, 0.00289343, -2.55565e-08, 1.91684e-08,
      3.01647e-09, 6.84051e-09, 2.79592e-10, -1.35767e-09, 0.00420618, 0.092451,
      -0.0239374, 0.0036276, 0.0113951, 0.0471896, 0.0465325, -0.00398807,
      -0.00484251, 0.0051614, -0.0031325, -0.0112447, -0.010968, -0.000896935,
      -0.00488789, 0.000951266, 0.0239709, -0.00630323, 0.000419006, 0.0201651,
      -0.00573095, -0.00118124, -0.0269275, 0.00700955, -0.00345217, 0.00356488,
      0.0134361, 0.013156, 9.58532e-05, 0.00315613, -2.58268e-05, -0.00124098,
      0.000352706, -1.80679e-06, -8.71959e-05, 2.47799e-05, -0.0150073,
      0.0114186, 0.0443012, 0.0180327, -0.026287, 0.0734351, -0.0643994,
      0.0257628, 0.0132084, -0.0123339, 0.0092297, -0.0148779, 0.0122493,
      -0.00246837, -0.0125735, -0.00375172, 0.00294872, 0.0112899, 0.00648758,
      -0.0055755, -0.0191436, 0.00433063, -0.00332619, -0.0128321, 0.0111166,
      -0.00969272, 0.0189659, -0.0160119, 0.00458659, 0.0107104, -0.000399315,
      0.000343129, 0.00117813, -2.80525e-05, 2.41086e-05, 8.2778e-05,
      -0.0450479, -0.00034974, -0.015063, 0.0655838, 0.0115163, -0.022358,
      0.0220978, 0.0583236, 0.0513097, -0.0119156, 0.00490159, 0.0116429,
      -0.0132479, -0.0146227, -0.00863565, -0.0118978, -0.000282044,
      -0.00400272, 0.0199347, 0.00139057, 0.00635067, 0.0131991, 0.000163177,
      0.00441453, 0.0159091, -0.00241207, -0.0110696, 0.0123057, 0.0171503,
      0.0119626, -0.00122682, -8.55716e-05, -0.00039083, -8.62007e-05,
      -6.0128e-06, -2.746e-05, -0.0304688, -0.954049, 0.259333, 0.0971056,
      0.255313, 0.950672, 0.983887, -0.0547431, -0.0857965, 0.0170489,
      -0.0104387, -0.036743, -0.0381557, -0.00278036, -0.0169143, -0.00177889,
      -0.04041, 0.0106552, -2.23782e-07, 2.40738e-07, 1.03401e-07, -0.000182535,
      -0.00516415, 0.00138942, 0.00125201, -0.00139237, -0.00501195,
      -0.00519809, -0.000154171, -0.00125602, 4.03664e-08, -6.04796e-08,
      -4.6768e-08, -2.38233e-09, 2.31605e-09, 1.35922e-09;

  A.BSESinglets().eigenvectors() = spsi_ref;

  // reference energy
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;

  A.BSESinglets().eigenvalues() = se_ref;
  A.CalcCoupledTransition_Dipoles();

  BOOST_CHECK_EQUAL(rho_f->NeedsInitialState(), true);

  rho_f->UpdateHist(A, QMState("s1"));

  std::vector<votca::Index> ref = {0};
  std::vector<votca::Index> results =
      rho_f->CalcIndeces(A, QMStateType::Singlet);
  for (votca::Index i = 0; i < votca::Index(ref.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i], results[i]);
  }

  BOOST_CHECK_EQUAL(results.size(), ref.size());

  rho_f->UpdateHist(A, QMState("s2"));
  std::vector<votca::Index> results2 =
      rho_f->CalcIndeces(A, QMStateType::Singlet);
  std::vector<votca::Index> ref2 = {1};
  for (votca::Index i = 0; i < votca::Index(ref2.size()); i++) {
    BOOST_CHECK_EQUAL(ref2[i], results2[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
