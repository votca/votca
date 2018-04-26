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

#define BOOST_TEST_MODULE eris_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/ERIs.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(eris_test)

BOOST_AUTO_TEST_CASE(fourcenter_direct){
  
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
 
 Eigen::MatrixXd dmat=Eigen::MatrixXd::Zero(17,17);       
 dmat<<0.00157507,0.0337454,4.48905e-16,-5.93152e-16,7.87133e-17,0.030876,2.51254e-16,-1.49094e-16,5.77899e-17,0.00415998,-0.00445632,0.00415998,-0.00445632,0.00415998,-0.00445632,0.00415998,-0.00445632,
0.0337454,0.722983,2.66427e-15,-4.44783e-15,3.45846e-16,0.661507,4.39854e-15,-2.02475e-15,1.04832e-15,0.0891262,-0.095475,0.0891262,-0.095475,0.0891262,-0.095475,0.0891262,-0.095475,
4.48905e-16,2.66427e-15,1.52199,2.88658e-15,2.09034e-15,-7.94212e-15,0.215492,2.8727e-15,-1.40513e-15,0.141933,-0.0402359,0.141933,-0.0402359,-0.141933,0.0402359,-0.141933,0.0402359,
-5.93152e-16,-4.44783e-15,2.88658e-15,1.52199,-2.31759e-15,9.21105e-15,-2.22045e-15,0.215492,1.6263e-15,0.141933,-0.0402359,-0.141933,0.0402359,-0.141933,0.0402359,0.141933,-0.0402359,
7.87133e-17,3.45846e-16,2.09034e-15,-2.31759e-15,1.52199,2.98902e-15,-2.04958e-15,4.79738e-15,0.215492,0.141933,-0.0402359,-0.141933,0.0402359,0.141933,-0.0402359,-0.141933,0.0402359,
0.030876,0.661507,-7.94212e-15,9.21105e-15,2.98902e-15,0.605259,2.55488e-15,2.7779e-17,1.33759e-15,0.0815477,-0.0873567,0.0815477,-0.0873567,0.0815477,-0.0873567,0.0815477,-0.0873567,
2.51254e-16,4.39854e-15,0.215492,-2.22045e-15,-2.04958e-15,2.55488e-15,0.0305108,3.29597e-17,-5.29036e-16,0.0200958,-0.00569686,0.0200958,-0.00569686,-0.0200958,0.00569686,-0.0200958,0.00569686,
-1.49094e-16,-2.02475e-15,2.8727e-15,0.215492,4.79738e-15,2.7779e-17,3.29597e-17,0.0305108,9.55941e-16,0.0200958,-0.00569686,-0.0200958,0.00569686,-0.0200958,0.00569686,0.0200958,-0.00569686,
5.77899e-17,1.04832e-15,-1.40513e-15,1.6263e-15,0.215492,1.33759e-15,-5.29036e-16,9.55941e-16,0.0305108,0.0200958,-0.00569686,-0.0200958,0.00569686,0.0200958,-0.00569686,-0.0200958,0.00569686,
0.00415998,0.0891262,0.141933,0.141933,0.141933,0.0815477,0.0200958,0.0200958,0.0200958,0.0506951,-0.0230264,-0.00224894,-0.00801753,-0.00224894,-0.00801753,-0.00224894,-0.00801753,
-0.00445632,-0.095475,-0.0402359,-0.0402359,-0.0402359,-0.0873567,-0.00569686,-0.00569686,-0.00569686,-0.0230264,0.0157992,-0.00801753,0.0115445,-0.00801753,0.0115445,-0.00801753,0.0115445,
0.00415998,0.0891262,0.141933,-0.141933,-0.141933,0.0815477,0.0200958,-0.0200958,-0.0200958,-0.00224894,-0.00801753,0.0506951,-0.0230264,-0.00224894,-0.00801753,-0.00224894,-0.00801753,
-0.00445632,-0.095475,-0.0402359,0.0402359,0.0402359,-0.0873567,-0.00569686,0.00569686,0.00569686,-0.00801753,0.0115445,-0.0230264,0.0157992,-0.00801753,0.0115445,-0.00801753,0.0115445,
0.00415998,0.0891262,-0.141933,-0.141933,0.141933,0.0815477,-0.0200958,-0.0200958,0.0200958,-0.00224894,-0.00801753,-0.00224894,-0.00801753,0.0506951,-0.0230264,-0.00224894,-0.00801753,
-0.00445632,-0.095475,0.0402359,0.0402359,-0.0402359,-0.0873567,0.00569686,0.00569686,-0.00569686,-0.00801753,0.0115445,-0.00801753,0.0115445,-0.0230264,0.0157992,-0.00801753,0.0115445,
0.00415998,0.0891262,-0.141933,0.141933,-0.141933,0.0815477,-0.0200958,0.0200958,-0.0200958,-0.00224894,-0.00801753,-0.00224894,-0.00801753,-0.00224894,-0.00801753,0.0506951,-0.0230264,
-0.00445632,-0.095475,0.0402359,-0.0402359,0.0402359,-0.0873567,0.00569686,-0.00569686,0.00569686,-0.00801753,0.0115445,-0.00801753,0.0115445,-0.00801753,0.0115445,-0.0230264,0.0157992;


ERIs eris1;
ERIs eris2;
eris1.Initialize_4c_screening(aobasis,1e-10);
eris2.Initialize_4c_small_molecule(aobasis);
eris1.CalculateERIs_4c_direct(aobasis,dmat);
eris2.CalculateERIs_4c_small_molecule(dmat);

bool check_eris=eris1.getERIs().isApprox(eris2.getERIs(),0.001);
if(!check_eris){
  std::cout<<eris1.getERIs()<<std::endl;
  std::cout<<eris2.getERIs()<<std::endl;
}
BOOST_CHECK_EQUAL(check_eris, 1);

}
BOOST_AUTO_TEST_SUITE_END()
