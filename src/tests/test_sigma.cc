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

#define BOOST_TEST_MODULE sigma_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/sigma.h>


using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full){
  
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
  rpa.configure(4,0,16);
  
  PPM ppm;
  Eigen::VectorXd screen_r=Eigen::VectorXd::Zero(1);
  screen_r(0)=ppm.getScreening_r();
  Eigen::VectorXd screen_i=Eigen::VectorXd::Zero(1);
  screen_i(0)=ppm.getScreening_i();
  rpa.setScreening(screen_r,screen_i);
  rpa.calculate_epsilon(es.eigenvalues(),Mmn);
  ppm.PPM_construct_parameters(rpa);
   Mmn.MultiplyLeftWithAuxMatrix(ppm.getPpm_phi_T());
  
  votca::ctp::Logger _log;
  Sigma sigma=Sigma(&_log);
  sigma.configure(4,0,16,5,0.001);
  Eigen::MatrixXd Vxc=Eigen::MatrixXd::Zero(17,17);
  sigma.setDFTdata(0.0,&Vxc,&es.eigenvalues());
  sigma.setGWAEnergies(es.eigenvalues());
  sigma.CalcdiagElements(Mmn,ppm);
  
  Eigen::VectorXd qp_energy_ref=Eigen::VectorXd::Zero(17);
  qp_energy_ref<< 0.164552, -0.21215, -0.21215,-0.276217, -21.5846, -14.6766, -14.6766, -17.8938, -31.0487, -73.2587, -73.2587, -205.47, -500.169, -255.833, -255.833, -383.092, -74.8316;
  bool qp_check=qp_energy_ref.isApprox(sigma.getGWAEnergies(),0.0001);
   if(!qp_check){
      cout<<"qp_energy"<<endl;
      cout<<sigma.getGWAEnergies()<<endl;
      cout<<"qp_energy_ref"<<endl;
      cout<<qp_energy_ref<<endl;
  }
  
  sigma.CalcOffDiagElements(Mmn,ppm);
  
  
  Eigen::MatrixXd Hqp=sigma.SetupFullQPHamiltonian();
  
  Eigen::MatrixXd Hqp_ref=Eigen::MatrixXd::Zero(17,17);
  Hqp_ref<<0.164552, 1.52509e-15, -1.64356e-15, -0.169138, -1.98272, -1.48687e-15, 4.9115e-15, -2.14042, -3.34583, 1.03505e-13, 8.97787e-14, -4.35298, 6.92053, -2.63768e-14, -6.80716e-15, -4.65247, 0.646687,
1.52509e-15, -0.21215, -5.06797e-15, 3.83751e-17, 2.34269e-15, 0.0118006, 0.78131, -2.07298e-15, -6.39938e-15, 3.43378, 1.42911, -9.49754e-16, -3.12305e-14, 1.91991, 0.150441, 8.21926e-15, -4.71319e-17,
-1.64356e-15, -5.06797e-15, -0.21215, 3.21142e-16, 7.91623e-16, -0.78131, 0.0118006, -9.93297e-16, 1.02404e-14, 1.42911, -3.43378, 4.92234e-14, -2.81844e-14, -0.150441, 1.91991, 4.12558e-15, -2.11895e-15,
-0.169138, 3.83751e-17, 3.21142e-16, -0.276217, -0.741458, 6.3439e-15, 4.75796e-15, 0.789749, 1.21449, -7.24139e-14, -3.88425e-14, 3.37089, 3.38614, 3.73137e-14, 1.23966e-14, 1.84539, -0.279895,
-1.98272, 2.34269e-15, 7.91623e-16, -0.741458, -21.5846, 1.75855e-14, 8.16739e-15, 1.88537, 10.2697, 4.09592e-13, 3.64945e-13, -49.8155, 94.9642, -2.17162e-13, 9.73485e-14, -4.22806, -2.33045,
-1.48687e-15, 0.0118006, -0.78131, 6.3439e-15, 1.75855e-14, -14.6766, 9.30428e-14, -7.51632e-15, -1.10574e-14, 2.97216, -7.45682, 6.46679e-14, -1.30112e-13, -3.15747, 33.7435, -5.73405e-14, 1.99509e-15,
4.9115e-15, 0.78131, 0.0118006, 4.75796e-15, 8.16739e-15, 9.30428e-14, -14.6766, 5.23429e-16, -1.75846e-14, -7.45682, -2.97216, 5.34149e-14, -1.19073e-13, -33.7435, -3.15747, 2.77312e-14, 6.53785e-15,
-2.14042, -2.07298e-15, -9.93297e-16, 0.789749, 1.88537, -7.51632e-15, 5.23429e-16, -17.8938, -0.955305, 1.62349e-13, 1.44457e-13, 0.0787898, 8.36091, -4.11584e-14, -4.89325e-14, -59.8084, -0.787966,
-3.34583, -6.39938e-15, 1.02404e-14, 1.21449, 10.2697, -1.10574e-14, -1.75846e-14, -0.955305, -31.0487, -9.22735e-13, -7.94657e-13, 6.30548, -5.30572, 5.42111e-14, 1.42055e-14, -24.1191, 2.4575,
1.03505e-13, 3.43378, 1.42911, -7.24139e-14, 4.09592e-13, 2.97216, -7.45682, 1.62349e-13, -9.22735e-13, -73.2587, -5.14588e-14, 9.61511e-13, -1.10547e-12, -53.3248, -27.2608, 1.72667e-13, 5.53087e-14,
8.97787e-14, 1.42911, -3.43378, -3.88425e-14, 3.64945e-13, -7.45682, -2.97216, 1.44457e-13, -7.94657e-13, -5.14588e-14, -73.2587, 9.05028e-13, -1.12346e-12, -27.2608, 53.3248, 1.32466e-13, 6.20531e-14,
-4.35298, -9.49754e-16, 4.92234e-14, 3.37089, -49.8155, 6.46679e-14, 5.34149e-14, 0.0787898, 6.30548, 9.61511e-13, 9.05028e-13, -205.47, 219.455, -3.17492e-13, 1.17636e-14, -75.4131, 14.0117,
6.92053, -3.12305e-14, -2.81844e-14, 3.38614, 94.9642, -1.30112e-13, -1.19073e-13, 8.36091, -5.30572, -1.10547e-12, -1.12346e-12, 219.455, -500.169, 4.44397e-13, -6.76416e-14, 92.9191, -19.1719,
-2.63768e-14, 1.91991, -0.150441, 3.73137e-14, -2.17162e-13, -3.15747, -33.7435, -4.11584e-14, 5.42111e-14, -53.3248, -27.2608, -3.17492e-13, 4.44397e-13, -255.833, 2.70814e-14, -1.25318e-13, 1.9266e-14,
-6.80716e-15, 0.150441, 1.91991, 1.23966e-14, 9.73485e-14, 33.7435, -3.15747, -4.89325e-14, 1.42055e-14, -27.2608, 53.3248, 1.17636e-14, -6.76416e-14, 2.70814e-14, -255.833, 4.67593e-14, -2.66148e-14,
-4.65247, 8.21926e-15, 4.12558e-15, 1.84539, -4.22806, -5.73405e-14, 2.77312e-14, -59.8084, -24.1191, 1.72667e-13, 1.32466e-13, -75.4131, 92.9191, -1.25318e-13, 4.67593e-14, -383.092, 12.917,
0.646687, -4.71319e-17, -2.11895e-15, -0.279895, -2.33045, 1.99509e-15, 6.53785e-15, -0.787966, 2.4575, 5.53087e-14, 6.20531e-14, 14.0117, -19.1719, 1.9266e-14, -2.66148e-14, 12.917, -74.8316;

  bool hqp_check=Hqp_ref.isApprox(Hqp,0.0001);
  if(!hqp_check){
      cout<<"hqp"<<endl;
      cout<<Hqp<<endl;
      cout<<"hqp_ref"<<endl;
      cout<<Hqp_ref<<endl;
  }
  
   BOOST_CHECK_EQUAL(hqp_check, true);
 BOOST_CHECK_EQUAL(qp_check, true);


}
        
BOOST_AUTO_TEST_SUITE_END()
