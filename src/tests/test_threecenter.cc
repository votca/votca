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

#define BOOST_TEST_MODULE threecenter_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/ERIs.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(threecenter_test)

BOOST_AUTO_TEST_CASE(threecenter_dft) {
  
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
  TCMatrix_dft threec;
  threec.Fill(aobasis,aobasis);
  
  
  
  
  Eigen::MatrixXd Res0=Eigen::MatrixXd::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
  threec.getDatamatrix(0).AddtoEigenMatrix(Res0);

  
  Eigen::MatrixXd Res4=Eigen::MatrixXd::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
  threec.getDatamatrix(4).AddtoEigenMatrix(Res4);

  
  Eigen::MatrixXd Ref0=Eigen::MatrixXd::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
  Ref0<<2.29374 ,0.304177 ,0 ,0 ,0 ,0.337068 ,0 ,0 ,0 ,0.033705 ,0.150367 ,0.033705 ,0.150367 ,0.033705 ,0.150367 ,0.033705 ,0.150367,
0.304177 ,0.921384 ,0 ,0 ,0 ,0.632988 ,0 ,0 ,0 ,0.138532 ,0.319933 ,0.138532 ,0.319933 ,0.138532 ,0.319933 ,0.138532 ,0.319933,
0 ,0 ,0.967997 ,0 ,0 ,0 ,0.366964 ,0 ,0 ,0.109996 ,0.0905438 ,0.109996 ,0.0905438 ,-0.109996 ,-0.0905438 ,-0.109996 ,-0.0905438,
0 ,0 ,0 ,0.967997 ,0 ,0 ,0 ,0.366964 ,0 ,0.109996 ,0.0905438 ,-0.109996 ,-0.0905438 ,-0.109996 ,-0.0905438 ,0.109996 ,0.0905438,
0 ,0 ,0 ,0 ,0.967997 ,0 ,0 ,0 ,0.366964 ,0.109996 ,0.0905438 ,-0.109996 ,-0.0905438 ,0.109996 ,0.0905438 ,-0.109996 ,-0.0905438,
0.337068 ,0.632988 ,0 ,0 ,0 ,0.616586 ,0 ,0 ,0 ,0.17704 ,0.362827 ,0.17704 ,0.362827 ,0.17704 ,0.362827 ,0.17704 ,0.362827,
0 ,0 ,0.366964 ,0 ,0 ,0 ,0.423084 ,0 ,0 ,0.130034 ,0.131792 ,0.130034 ,0.131792 ,-0.130034 ,-0.131792 ,-0.130034 ,-0.131792,
0 ,0 ,0 ,0.366964 ,0 ,0 ,0 ,0.423084 ,0 ,0.130034 ,0.131792 ,-0.130034 ,-0.131792 ,-0.130034 ,-0.131792 ,0.130034 ,0.131792,
0 ,0 ,0 ,0 ,0.366964 ,0 ,0 ,0 ,0.423084 ,0.130034 ,0.131792 ,-0.130034 ,-0.131792 ,0.130034 ,0.131792 ,-0.130034 ,-0.131792,
0.033705 ,0.138532 ,0.109996 ,0.109996 ,0.109996 ,0.17704 ,0.130034 ,0.130034 ,0.130034 ,0.436704 ,0.280477 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251,
0.150367 ,0.319933 ,0.0905438 ,0.0905438 ,0.0905438 ,0.362827 ,0.131792 ,0.131792 ,0.131792 ,0.280477 ,0.39915 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028,
0.033705 ,0.138532 ,0.109996 ,-0.109996 ,-0.109996 ,0.17704 ,0.130034 ,-0.130034 ,-0.130034 ,0.00555884 ,0.0631251 ,0.436704 ,0.280477 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251,
0.150367 ,0.319933 ,0.0905438 ,-0.0905438 ,-0.0905438 ,0.362827 ,0.131792 ,-0.131792 ,-0.131792 ,0.0631251 ,0.182028 ,0.280477 ,0.39915 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028,
0.033705 ,0.138532 ,-0.109996 ,-0.109996 ,0.109996 ,0.17704 ,-0.130034 ,-0.130034 ,0.130034 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251 ,0.436704 ,0.280477 ,0.00555884 ,0.0631251,
0.150367 ,0.319933 ,-0.0905438 ,-0.0905438 ,0.0905438 ,0.362827 ,-0.131792 ,-0.131792 ,0.131792 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028 ,0.280477 ,0.39915 ,0.0631251 ,0.182028,
0.033705 ,0.138532 ,-0.109996 ,0.109996 ,-0.109996 ,0.17704 ,-0.130034 ,0.130034 ,-0.130034 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251 ,0.00555884 ,0.0631251 ,0.436704 ,0.280477,
0.150367 ,0.319933 ,-0.0905438 ,0.0905438 ,-0.0905438 ,0.362827 ,-0.131792 ,0.131792 ,-0.131792 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028 ,0.0631251 ,0.182028 ,0.280477 ,0.39915;
  
  Eigen::MatrixXd Ref4=Eigen::MatrixXd::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
  Ref4<<0 ,0 ,0 ,0 ,0.21851 ,0 ,0 ,0 ,0.0305022 ,0.00646054 ,0.00674129 ,-0.00646054 ,-0.00674129 ,0.00646054 ,0.00674129 ,-0.00646054 ,-0.00674129,
0 ,0 ,0 ,0 ,0.896593 ,0 ,0 ,0 ,0.409428 ,0.125791 ,0.10102 ,-0.125791 ,-0.10102 ,0.125791 ,0.10102 ,-0.125791 ,-0.10102,
0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.0866627 ,0.0226036 ,-0.0866627 ,-0.0226036 ,-0.0866627 ,-0.0226036 ,0.0866627 ,0.0226036,
0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.0866627 ,0.0226036 ,0.0866627 ,0.0226036 ,-0.0866627 ,-0.0226036 ,-0.0866627 ,-0.0226036,
0.21851 ,0.896593 ,0 ,0 ,0 ,0.690122 ,0 ,0 ,0 ,0.17602 ,0.363256 ,0.17602 ,0.363256 ,0.17602 ,0.363256 ,0.17602 ,0.363256,
0 ,0 ,0 ,0 ,0.690122 ,0 ,0 ,0 ,0.56523 ,0.18263 ,0.159752 ,-0.18263 ,-0.159752 ,0.18263 ,0.159752 ,-0.18263 ,-0.159752,
0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.122829 ,0.0489282 ,-0.122829 ,-0.0489282 ,-0.122829 ,-0.0489282 ,0.122829 ,0.0489282,
0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0.122829 ,0.0489282 ,0.122829 ,0.0489282 ,-0.122829 ,-0.0489282 ,-0.122829 ,-0.0489282,
0.0305022 ,0.409428 ,0 ,0 ,0 ,0.56523 ,0 ,0 ,0 ,0.201612 ,0.37355 ,0.201612 ,0.37355 ,0.201612 ,0.37355 ,0.201612 ,0.37355,
0.00646054 ,0.125791 ,0.0866627 ,0.0866627 ,0.17602 ,0.18263 ,0.122829 ,0.122829 ,0.201612 ,0.540893 ,0.311822 ,2.46557e-19 ,0.0523107 ,0.00882417 ,0.0808349 ,2.44863e-19 ,0.0523107,
0.00674129 ,0.10102 ,0.0226036 ,0.0226036 ,0.363256 ,0.159752 ,0.0489282 ,0.0489282 ,0.37355 ,0.311822 ,0.304029 ,-0.0523107 ,3.0974e-18 ,0.0808349 ,0.15978 ,-0.0523107 ,3.0974e-18,
-0.00646054 ,-0.125791 ,-0.0866627 ,0.0866627 ,0.17602 ,-0.18263 ,-0.122829 ,0.122829 ,0.201612 ,2.46557e-19 ,-0.0523107 ,-0.540893 ,-0.311822 ,-2.44863e-19 ,-0.0523107 ,-0.00882417 ,-0.0808349,
-0.00674129 ,-0.10102 ,-0.0226036 ,0.0226036 ,0.363256 ,-0.159752 ,-0.0489282 ,0.0489282 ,0.37355 ,0.0523107 ,3.0974e-18 ,-0.311822 ,-0.304029 ,0.0523107 ,-3.0974e-18 ,-0.0808349 ,-0.15978,
0.00646054 ,0.125791 ,-0.0866627 ,-0.0866627 ,0.17602 ,0.18263 ,-0.122829 ,-0.122829 ,0.201612 ,0.00882417 ,0.0808349 ,-2.44863e-19 ,0.0523107 ,0.540893 ,0.311822 ,2.46557e-19 ,0.0523107,
0.00674129 ,0.10102 ,-0.0226036 ,-0.0226036 ,0.363256 ,0.159752 ,-0.0489282 ,-0.0489282 ,0.37355 ,0.0808349 ,0.15978 ,-0.0523107 ,-3.0974e-18 ,0.311822 ,0.304029 ,-0.0523107 ,3.0974e-18,
-0.00646054 ,-0.125791 ,0.0866627 ,-0.0866627 ,0.17602 ,-0.18263 ,0.122829 ,-0.122829 ,0.201612 ,2.44863e-19 ,-0.0523107 ,-0.00882417 ,-0.0808349 ,2.46557e-19 ,-0.0523107 ,-0.540893 ,-0.311822,
-0.00674129 ,-0.10102 ,0.0226036 ,-0.0226036 ,0.363256 ,-0.159752 ,0.0489282 ,-0.0489282 ,0.37355 ,0.0523107 ,3.0974e-18 ,-0.0808349 ,-0.15978 ,0.0523107 ,3.0974e-18 ,-0.311822 ,-0.304029;
  
  bool check_threec=(Res0.isApprox(Ref0,0.00001) && Res4.isApprox(Ref4,0.00001) );
BOOST_CHECK_EQUAL(check_threec, 1);

}

BOOST_AUTO_TEST_CASE(threecenter_gwbse){
  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(&basis,orbitals.QMAtoms());
  
  AOKinetic kinetic;
kinetic.Fill(aobasis);

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kinetic.Matrix());

TCMatrix_gwbse tc;
tc.Initialize(aobasis.AOBasisSize(),0,5,0,7);
tc.Fill(aobasis,aobasis,es.eigenvectors());
  
 
  
}
BOOST_AUTO_TEST_SUITE_END()
