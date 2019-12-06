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

#define BOOST_TEST_MODULE statetracker_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/statetracker.h>

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(statetracker_test)
BOOST_AUTO_TEST_CASE(osc) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orb;
  orb.QMAtoms().LoadFromFile("molecule.xyz");
  orb.setDFTbasisName("3-21G.xml");
  Logger log;
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);
  orb.setBSEindices(0, 16);
  orb.setTDAApprox(true);

  Eigen::MatrixXd& MOs = orb.MOs().eigenvectors();
  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  MOs = Eigen::MatrixXd::Zero(17, 17);
  MOs << -0.00761992, -4.69664e-13, 8.35009e-15, -1.15214e-14, -0.0156169,
      -2.23157e-12, 1.52916e-14, 2.10997e-15, 8.21478e-15, 3.18517e-15,
      2.89043e-13, -0.00949189, 1.95787e-12, 1.22168e-14, -2.63092e-15,
      -0.22227, 1.00844, 0.233602, -3.18103e-12, 4.05093e-14, -4.70943e-14,
      0.1578, 4.75897e-11, -1.87447e-13, -1.02418e-14, 6.44484e-14, -2.6602e-14,
      6.5906e-12, -0.281033, -6.67755e-12, 2.70339e-14, -9.78783e-14, -1.94373,
      -0.36629, -1.63678e-13, -0.22745, -0.054851, 0.30351, 3.78688e-11,
      -0.201627, -0.158318, -0.233561, -0.0509347, -0.650424, 0.452606,
      -5.88565e-11, 0.453936, -0.165715, -0.619056, 7.0149e-12, 2.395e-14,
      -4.51653e-14, -0.216509, 0.296975, -0.108582, 3.79159e-11, -0.199301,
      0.283114, -0.0198557, 0.584622, 0.275311, 0.461431, -5.93732e-11,
      0.453057, 0.619523, 0.166374, 7.13235e-12, 2.56811e-14, -9.0903e-14,
      -0.21966, -0.235919, -0.207249, 3.75979e-11, -0.199736, -0.122681,
      0.255585, -0.534902, 0.362837, 0.461224, -5.91028e-11, 0.453245,
      -0.453298, 0.453695, 7.01644e-12, 2.60987e-14, 0.480866, 1.8992e-11,
      -2.56795e-13, 4.14571e-13, 2.2709, 4.78615e-10, -2.39153e-12,
      -2.53852e-13, -2.15605e-13, -2.80359e-13, 7.00137e-12, 0.145171,
      -1.96136e-11, -2.24876e-13, -2.57294e-14, 4.04176, 0.193617, -1.64421e-12,
      -0.182159, -0.0439288, 0.243073, 1.80753e-10, -0.764779, -0.600505,
      -0.885907, 0.0862014, 1.10077, -0.765985, 6.65828e-11, -0.579266,
      0.211468, 0.789976, -1.41532e-11, -1.29659e-13, -1.64105e-12, -0.173397,
      0.23784, -0.0869607, 1.80537e-10, -0.755957, 1.07386, -0.0753135,
      -0.989408, -0.465933, -0.78092, 6.72256e-11, -0.578145, -0.790571,
      -0.212309, -1.42443e-11, -1.31306e-13, -1.63849e-12, -0.17592, -0.188941,
      -0.165981, 1.79403e-10, -0.757606, -0.465334, 0.969444, 0.905262,
      -0.61406, -0.78057, 6.69453e-11, -0.578385, 0.578453, -0.578959,
      -1.40917e-11, -1.31002e-13, 0.129798, -0.274485, 0.00256652, -0.00509635,
      -0.0118465, 0.141392, -0.000497905, -0.000510338, -0.000526798,
      -0.00532572, 0.596595, 0.65313, -0.964582, -0.000361559, -0.000717866,
      -0.195084, 0.0246232, 0.0541331, -0.255228, 0.00238646, -0.0047388,
      -0.88576, 1.68364, -0.00592888, -0.00607692, -9.5047e-05, -0.000960887,
      0.10764, -0.362701, 1.53456, 0.000575205, 0.00114206, -0.793844,
      -0.035336, 0.129798, 0.0863299, -0.0479412, 0.25617, -0.0118465,
      -0.0464689, 0.0750316, 0.110468, -0.0436647, -0.558989, -0.203909,
      0.65313, 0.320785, 0.235387, 0.878697, -0.195084, 0.0246232, 0.0541331,
      0.0802732, -0.0445777, 0.238198, -0.88576, -0.553335, 0.893449, 1.31541,
      -0.00787816, -0.100855, -0.0367902, -0.362701, -0.510338, -0.374479,
      -1.39792, -0.793844, -0.035336, 0.129798, 0.0927742, -0.197727, -0.166347,
      -0.0118465, -0.0473592, 0.0582544, -0.119815, -0.463559, 0.320126,
      -0.196433, 0.65313, 0.321765, 0.643254, -0.642737, -0.195084, 0.0246232,
      0.0541331, 0.0862654, -0.183855, -0.154677, -0.88576, -0.563936, 0.693672,
      -1.42672, -0.0836372, 0.0577585, -0.0354411, -0.362701, -0.511897,
      -1.02335, 1.02253, -0.793844, -0.035336, 0.129798, 0.0953806, 0.243102,
      -0.0847266, -0.0118465, -0.0475639, -0.132788, 0.00985812, 0.507751,
      0.244188, -0.196253, 0.65313, 0.322032, -0.87828, -0.235242, -0.195084,
      0.0246232, 0.0541331, 0.088689, 0.226046, -0.0787824, -0.88576, -0.566373,
      -1.58119, 0.117387, 0.0916104, 0.0440574, -0.0354087, -0.362701,
      -0.512321, 1.39726, 0.374248, -0.793844, -0.035336;

  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;
  orb.BSESinglets().eigenvalues() = se_ref;
  // reference coefficients
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
  orb.BSESinglets().eigenvectors() = spsi_ref;
  orb.CalcCoupledTransition_Dipoles();

  {
    std::ofstream options("statetracker.xml");
    options << "<statetracker>" << std::endl;
    options << "  <filters>oscillatorstrength</filters>" << std::endl;
    options << "  <oscillatorstrength>0.0019</oscillatorstrength>" << std::endl;
    options << "</statetracker>" << std::endl;
    options.close();
    votca::tools::Property prop;
    prop.LoadFromXML("statetracker.xml");
    StateTracker tracker;
    tracker.setLogger(&log);
    QMState s("s1");
    tracker.setInitialState(s);
    tracker.Initialize(prop.get("statetracker"));
    QMState newstate = tracker.CalcState(orb);
    BOOST_CHECK_EQUAL(newstate.Type().ToString(), "s");
    BOOST_CHECK_EQUAL(newstate.StateIdx(), 1);
  }

  {
    std::ofstream options("statetracker.xml");
    options << "<statetracker>" << std::endl;
    options << "  <filters>localisation</filters>" << std::endl;
    options << "  <localisation>" << std::endl;
    options << "  <fragment>0 1</fragment>" << std::endl;
    options << "  <threshold>0.5</threshold>" << std::endl;
    options << "  </localisation>" << std::endl;
    options << "</statetracker>" << std::endl;
    options.close();
    votca::tools::Property prop;
    prop.LoadFromXML("statetracker.xml");
    StateTracker tracker;
    tracker.setLogger(&log);
    QMState s("s1");
    tracker.setInitialState(s);
    tracker.Initialize(prop.get("statetracker"));
    QMState newstate = tracker.CalcState(orb);
    BOOST_CHECK_EQUAL(newstate.Type().ToString(), "s");
    BOOST_CHECK_EQUAL(newstate.StateIdx(), 0);
  }
}

BOOST_AUTO_TEST_CASE(readwrite_hdf5) {
  std::ofstream options("statetracker.xml");
  options << "<statetracker>" << std::endl;
  options << "  <filters>localisation</filters>" << std::endl;
  options << "  <localisation>" << std::endl;
  options << "  <fragment>0 1</fragment>" << std::endl;
  options << "  <threshold>0.1</threshold>" << std::endl;
  options << "  </localisation>" << std::endl;
  options << "</statetracker>" << std::endl;
  options.close();
  votca::tools::Property prop;
  prop.LoadFromXML("statetracker.xml");
  StateTracker tracker;
  QMState s("s1");
  tracker.setInitialState(s);
  tracker.Initialize(prop.get("statetracker"));

  {
    CheckpointFile f("statetracker_test.hdf5");
    CheckpointWriter w = f.getWriter();
    tracker.WriteToCpt(w);
  }
  StateTracker tracker2;
  CheckpointFile f("statetracker_test.hdf5");
  CheckpointReader r = f.getReader();
  tracker2.ReadFromCpt(r);

  BOOST_CHECK_EQUAL(tracker.InitialState().ToString(),
                    tracker2.InitialState().ToString());
}

BOOST_AUTO_TEST_SUITE_END()
