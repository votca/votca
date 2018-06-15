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

#define BOOST_TEST_MODULE ppm_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/ppm.h>


using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(ppm_test)

BOOST_AUTO_TEST_CASE(ppm_full){
  
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
  aobasis.AOBasisFill(&basis,orbitals.QMAtoms());
  
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
  rpa.configure(4,0,17-1);
  
  PPM ppm;
  Eigen::VectorXd screen_r=Eigen::VectorXd::Zero(1);
  screen_r(0)=ppm.getScreening_r();
  Eigen::VectorXd screen_i=Eigen::VectorXd::Zero(1);
  screen_i(0)=ppm.getScreening_i();
  rpa.setScreening(screen_r,screen_i);
  rpa.calculate_epsilon(es.eigenvalues(),Mmn);
  ppm.PPM_construct_parameters(rpa);
 
  
  Eigen::VectorXd ppm_freq=Eigen::VectorXd::Zero(17);
  ppm_freq<< 1.94598,0.923413,0.923413, 1.66455, 1.25616, 1.56824,0.939509,0.725025,0.725025, 0.69766, 0.69766,0.789006,0.655588,0.697099, 1.21575, 1.21575, 25.747;
  Eigen::VectorXd ppm_w=Eigen::VectorXd::Zero(17);
  ppm_w<<0.000145438,0.000350387,0.000350387,0.000393263,0.000511193, 0.00208969, 0.00455154, 0.0107557, 0.0107557, 0.0456673, 0.0456673, 0.089354, 0.118543, 0.788686, 0.984015, 0.984015, 0.994573;
 
  
  bool f_check =ppm_freq.isApprox(ppm.getPpm_freq(),0.0001);
  
  bool w_check =ppm_w.isApprox(ppm.getPpm_weight(),0.0001);
   
 
  
  if(!f_check){
      cout<<"ppm_freq"<<endl;
      cout<<ppm.getPpm_freq()<<endl;
      cout<<"ppm_freq_ref"<<endl;
      cout<<ppm_freq<<endl;
  }
 
  if(!w_check){
      cout<<"ppm_w"<<endl;
      cout<<ppm.getPpm_weight()<<endl;
      cout<<"ppm_w_ref"<<endl;
      cout<<ppm_w<<endl;
  }
 
 BOOST_CHECK_EQUAL(f_check, 1);
 BOOST_CHECK_EQUAL(w_check, 1);

}
        
BOOST_AUTO_TEST_SUITE_END()
