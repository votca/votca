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

#define BOOST_TEST_MODULE aobasis_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/convergenceacc.h>
#include <fstream>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(aobasis_test)

BOOST_AUTO_TEST_CASE(FillNormBasis_test) {   
std::ofstream basisfile("notnormalized.xml");
basisfile << "<basis name=\"def2-TZVP\">"<< std::endl;
basisfile << "  <element name=\"Al\">" << std::endl;
basisfile << "    <shell scale=\"1.0\" type=\"D\">"<< std::endl;
basisfile << "      <constant decay=\"1.570000e+00\">"<< std::endl;
basisfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"<<std::endl;
basisfile << "      </constant>"<<std::endl;
basisfile << "      <constant decay=\"3.330000e-01\">"<<std::endl;
basisfile << "        <contractions factor=\"1.000000e+00\" type=\"D\"/>"<<std::endl;
basisfile << "      </constant>"<<std::endl;
basisfile << "    </shell> "<<std::endl;
basisfile << "  </element> "<<std::endl;
basisfile << "</basis> "<<std::endl;
basisfile.close();
   
std::ofstream xyzfile("Al.xyz");
xyzfile << " 1" << std::endl;
xyzfile << " Al" << std::endl;
xyzfile << " Al            .000000     .000000     .000000" << std::endl;
xyzfile.close();

Orbitals orbitals;
orbitals.LoadFromXYZ("Al.xyz");
BasisSet basis;
basis.LoadBasisSet("notnormalized.xml");
AOBasis aobasis;
aobasis.AOBasisFill(basis,orbitals.QMAtoms());

const AOShell* shell=aobasis.getShell(0);
AOShell::GaussianIterator it;
std::vector<double> ref_results={0.1831079647,0.9155398233};
int i=0;
bool check_norm=true;
for(it=shell->firstGaussian();it<shell->lastGaussian();++it){
  if(std::abs(ref_results[i]-it->getContraction()[2])>1e-7){
   check_norm=false;
   break;
  }
  i++;
}


i=0;
if(!check_norm){
  for(it=shell->firstGaussian();it<shell->lastGaussian();++it){
  std::cout<<"Ref:"<<ref_results[i]<<" result:"<<it->getContraction()[2]<<std::endl;
   i++;
  }
 
}
BOOST_CHECK_EQUAL(check_norm, 1);

}


BOOST_AUTO_TEST_CASE(ReorderMos_test) {
  
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
  basisfile << "    <shell scale=\"1.0\" type=\"D\">"<< std::endl;
  basisfile << "      <constant decay=\"1.570000e+00\">"<< std::endl;
  basisfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"<<std::endl;
  basisfile << "      </constant>"<<std::endl;
  basisfile << "    </shell> "<<std::endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();
  
  Orbitals orbitals;
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(basis,orbitals.QMAtoms());
  AOOverlap overlap;
  overlap.Fill(aobasis);
    
  AOKinetic kinetic;
  kinetic.Fill(aobasis);
  AOESP esp;
  esp.Fillnucpotential(aobasis,orbitals.QMAtoms());
  Eigen::MatrixXd H=kinetic.Matrix()+esp.getNuclearpotential();
  
  double levelshift=0.1;
  ConvergenceAcc d;
  Orbitals orb;
  int occlevels=5;
  d.Configure(ConvergenceAcc::closed,false,false,10,false,0,0,levelshift,0,occlevels,0);
  d.setOverlap(&overlap.Matrix(),1e-8);
  d.SolveFockmatrix(orb.MOEnergies(),orb.MOCoefficients(),H);
  
  Eigen::MatrixXd ref=orb.MOCoefficients();
  aobasis.ReorderMOs(orb.MOCoefficients(),"xtp","gaussian");
  aobasis.ReorderMOs(orb.MOCoefficients(),"gaussian","xtp");
  
  
  Eigen::MatrixXd overlap_b=overlap.Matrix();
  aobasis.ReorderMatrix(overlap_b,"xtp","nwchem");
  aobasis.ReorderMatrix(overlap_b,"nwchem","xtp");
  bool check_reorder_overlap=overlap_b.isApprox(overlap.Matrix(),1e-7);
  BOOST_CHECK_EQUAL(check_reorder_overlap, 1);
  
  bool check_reorder_gaus=ref.isApprox(orb.MOCoefficients(),1e-7);
  if(!check_reorder_gaus){
    std::cout<<"ref"<<std::endl;
    std::cout<<ref<<std::endl;
    std::cout<<"reordered"<<std::endl;
    std::cout<<orb.MOCoefficients()<<std::endl;
  }
  
  BOOST_CHECK_EQUAL(check_reorder_gaus, 1);
  
  aobasis.ReorderMOs(orb.MOCoefficients(),"xtp","nwchem");
  aobasis.ReorderMOs(orb.MOCoefficients(),"nwchem","xtp");
  
   bool check_reorder_nwchem=ref.isApprox(orb.MOCoefficients(),1e-7);
  
  
  BOOST_CHECK_EQUAL(check_reorder_nwchem, 1);
}   


BOOST_AUTO_TEST_SUITE_END()
