/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/test/unit_test.hpp>
#include <votca/xtp/rpa.h>


using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_full){
  
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
  basisfile <<"<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>" << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();
  
  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  
  AOBasis aobasis;
  aobasis.AOBasisFill(basis,orbitals.QMAtoms());
  
  Orbitals orb;
  orb.setBasisSetSize(17);
  orb.setNumberOfLevels(4,13);

AOKinetic kinetic;
kinetic.Fill(aobasis);

AOESP esp;
esp.Fillnucpotential(aobasis,orbitals.QMAtoms());

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kinetic.Matrix()+esp.Matrix());

TCMatrix_gwbse Mmn;
Mmn.Initialize(aobasis.AOBasisSize(),0,16,0,16);
Mmn.Fill(aobasis,aobasis,es.eigenvectors());

  RPA rpa;
  rpa.configure(4,0,16);
  Eigen::VectorXd screen_r=Eigen::VectorXd::Zero(1);
  screen_r(0)=0;
  Eigen::VectorXd screen_i=Eigen::VectorXd::Zero(1);
  screen_i(0)=0.5;
  rpa.setScreening(screen_r,screen_i);
  rpa.calculate_epsilon(es.eigenvalues(),Mmn);
  
  Eigen::MatrixXd i_ref=Eigen::MatrixXd::Zero(17,17);
  i_ref<<1.20514,1.19201,0.13208,-0.13208,0.13208,2.59243,0.360147,-0.360147,0.360147,0.94804,2.53836,0.94804,2.53836,0.94804,2.53836,0.624023,1.97651,
1.19201,7.9361,0.769578,-0.769578,0.769578,15.0954,2.09638,-2.09638,2.09638,5.52074,14.7852,5.52074,14.7852,5.52074,14.7852,3.63538,11.516,
0.13208,0.769578,1.25496,-0.0357886,0.0357886,1.67041,0.768113,-0.0590628,0.0590628,0.751575,1.87312,0.751575,1.87312,0.358466,1.22231,0.370614,1.2096,
-0.13208,-0.769578,-0.0357886,1.25496,-0.0357886,-1.67041,-0.0590628,0.768113,-0.0590628,-0.358466,-1.22231,-0.751575,-1.87312,-0.751575,-1.87312,-0.370614,-1.2096,
0.13208,0.769578,0.0357886,-0.0357886,1.25496,1.67041,0.0590628,-0.0590628,0.768113,0.751575,1.87312,0.358466,1.22231,0.751575,1.87312,0.370614,1.2096,
2.59243,15.0954,1.67041,-1.67041,1.67041,33.9617,4.5658,-4.5658,4.5658,12.0692,32.3421,12.0692,32.3421,12.0692,32.3421,7.95011,25.2075,
0.360147,2.09638,0.768113,-0.0590628,0.0590628,4.5658,3.40668,-0.0454868,0.0454868,2.15163,5.288,2.15163,5.288,0.825447,3.08396,0.984063,3.25425,
-0.360147,-2.09638,-0.0590628,0.768113,-0.0590628,-4.5658,-0.0454868,3.40668,-0.0454868,-0.825447,-3.08396,-2.15163,-5.288,-2.15163,-5.288,-0.984063,-3.25425,
0.360147,2.09638,0.0590628,-0.0590628,0.768113,4.5658,0.0454868,-0.0454868,3.40668,2.15163,5.288,0.825447,3.08396,2.15163,5.288,0.984063,3.25425,
0.94804,5.52074,0.751575,-0.358466,0.751575,12.0692,2.15163,-0.825447,2.15163,5.94322,12.7086,4.17617,11.4491,4.17617,11.4491,2.88434,9.18305,
2.53836,14.7852,1.87312,-1.22231,1.87312,32.3421,5.288,-3.08396,5.288,12.7086,34.1927,11.4491,31.106,11.4491,31.106,7.75321,24.6616,
0.94804,5.52074,0.751575,-0.751575,0.358466,12.0692,2.15163,-2.15163,0.825447,4.17617,11.4491,5.94322,12.7086,4.17617,11.4491,2.88434,9.18305,
2.53836,14.7852,1.87312,-1.87312,1.22231,32.3421,5.288,-5.288,3.08396,11.4491,31.106,12.7086,34.1927,11.4491,31.106,7.75321,24.6616,
0.94804,5.52074,0.358466,-0.751575,0.751575,12.0692,0.825447,-2.15163,2.15163,4.17617,11.4491,4.17617,11.4491,5.94322,12.7086,2.88434,9.18305,
2.53836,14.7852,1.22231,-1.87312,1.87312,32.3421,3.08396,-5.288,5.288,11.4491,31.106,11.4491,31.106,12.7086,34.1927,7.75321,24.6616,
0.624023,3.63538,0.370614,-0.370614,0.370614,7.95011,0.984063,-0.984063,0.984063,2.88434,7.75321,2.88434,7.75321,2.88434,7.75321,3.00374,6.24909,
1.97651,11.516,1.2096,-1.2096,1.2096,25.2075,3.25425,-3.25425,3.25425,9.18305,24.6616,9.18305,24.6616,9.18305,24.6616,6.24909,20.6288;
  
  bool i_check =i_ref.isApprox(rpa.GetEpsilon_i()[0],0.0001);
  
   if(!i_check){
  cout<<"Epsilon_i"<<endl;
   cout<<rpa.GetEpsilon_i()[0]<<endl;
   cout<<"Epsilon_i_ref"<<endl;
   cout<<i_ref<<endl;
 }
 BOOST_CHECK_EQUAL(i_check , 1);
   
 Eigen::MatrixXd r_ref=Eigen::MatrixXd::Zero(17,17);
 r_ref<<1.21743,1.26312,0.1442,-0.1442,0.1442,2.74918,0.395553,-0.395553,0.395553,1.00786,2.69762,1.00786,2.69762,1.00786,2.69762,0.655853,2.08334,
1.26312,8.34802,0.838912,-0.838912,0.838912,16.0064,2.29882,-2.29882,2.29882,5.86828,15.7118,5.86828,15.7118,5.86828,15.7118,3.8217,12.1415,
0.1442,0.838912,2.89939,0.664883,-0.664883,1.82223,6.28539,2.35111,-2.35111,2.10307,4.41399,2.10307,4.41399,-2.09496,-3.22344,0.331895,1.14524,
-0.1442,-0.838912,0.664883,2.89939,0.664883,-1.82223,2.35111,6.28539,2.35111,2.09496,3.22344,-2.10307,-4.41399,-2.10307,-4.41399,-0.331895,-1.14524,
0.1442,0.838912,-0.664883,0.664883,2.89939,1.82223,-2.35111,2.35111,6.28539,2.10307,4.41399,-2.09496,-3.22344,2.10307,4.41399,0.331895,1.14524,
2.74918,16.0064,1.82223,-1.82223,1.82223,36.007,5.01152,-5.01152,5.01152,12.8525,34.4397,12.8525,34.4397,12.8525,34.4397,8.37311,26.6357,
0.395553,2.29882,6.28539,2.35111,-2.35111,5.01152,21.9541,8.21795,-8.21795,6.68798,13.7921,6.68798,13.7921,-7.50907,-12.0483,0.865854,3.06061,
-0.395553,-2.29882,2.35111,6.28539,2.35111,-5.01152,8.21795,21.9541,8.21795,7.50907,12.0483,-6.68798,-13.7921,-6.68798,-13.7921,-0.865854,-3.06061,
0.395553,2.29882,-2.35111,2.35111,6.28539,5.01152,-8.21795,8.21795,21.9541,6.68798,13.7921,-7.50907,-12.0483,6.68798,13.7921,0.865854,3.06061,
1.00786,5.86828,2.10307,2.09496,2.10307,12.8525,6.68798,7.50907,6.68798,10.3766,21.1055,2.42655,8.49486,2.42655,8.49486,3.01594,9.66101,
2.69762,15.7118,4.41399,3.22344,4.41399,34.4397,13.7921,12.0483,13.7921,21.1055,50.3046,8.49486,26.3736,8.49486,26.3736,8.12144,25.9812,
1.00786,5.86828,2.10307,-2.10307,-2.09496,12.8525,6.68798,-6.68798,-7.50907,2.42655,8.49486,10.3766,21.1055,2.42655,8.49486,3.01594,9.66101,
2.69762,15.7118,4.41399,-4.41399,-3.22344,34.4397,13.7921,-13.7921,-12.0483,8.49486,26.3736,21.1055,50.3046,8.49486,26.3736,8.12144,25.9812,
1.00786,5.86828,-2.09496,-2.10307,2.10307,12.8525,-7.50907,-6.68798,6.68798,2.42655,8.49486,2.42655,8.49486,10.3766,21.1055,3.01594,9.66101,
2.69762,15.7118,-3.22344,-4.41399,4.41399,34.4397,-12.0483,-13.7921,13.7921,8.49486,26.3736,8.49486,26.3736,21.1055,50.3046,8.12144,25.9812,
0.655853,3.8217,0.331895,-0.331895,0.331895,8.37311,0.865854,-0.865854,0.865854,3.01594,8.12144,3.01594,8.12144,3.01594,8.12144,3.17835,6.74615,
2.08334,12.1415,1.14524,-1.14524,1.14524,26.6357,3.06061,-3.06061,3.06061,9.66101,25.9812,9.66101,25.9812,9.66101,25.9812,6.74615,22.1044;
 bool r_check =r_ref.isApprox(rpa.GetEpsilon_r()[0],0.0001);
 
 
 if(!r_check){
    cout<<"Epsilon_r"<<endl;
   cout<<rpa.GetEpsilon_r()[0]<<endl;
    cout<<"Epsilon_r_ref"<<endl;
   cout<<r_ref<<endl;
 }

 BOOST_CHECK_EQUAL(r_check , 1);

}
        
BOOST_AUTO_TEST_SUITE_END()
