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

#define BOOST_TEST_MODULE polararsite_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/polarsite.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(polararsite_test)

BOOST_AUTO_TEST_CASE(constructors_test) { PolarSite ps(1, "ps1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  PolarSite ps(1,"ps2");
  BOOST_CHECK_EQUAL(ps.getId(),1);
  BOOST_CHECK_EQUAL(ps.getName(),"ps2");
}

BOOST_AUTO_TEST_CASE(multipole_test) {
  PolarSite ps(1,"ps2");
  Eigen::VectorXd multipoles=Eigen::VectorXd::Zero(9);
  multipoles<<1,2,3,4,8,7,2,3.3,-0.5;
  ps.setMultipoles(multipoles);
  bool check_mpoles=multipoles.isApprox(ps.getMultipoles(),0.0001);
   BOOST_CHECK_EQUAL(check_mpoles,true);
   
   bool check_rank=(ps.getRank()==2);
   BOOST_CHECK_EQUAL(check_rank,true);
  
}

BOOST_AUTO_TEST_CASE(translate_test) {
  PolarSite ps(1,"ps2");
  Eigen::Vector3d shift;
  shift<<0,0,5;
  ps.Translate(shift);
  BOOST_CHECK_EQUAL(shift.isApprox(ps.getPos(),1e-5),true);
}


BOOST_AUTO_TEST_CASE(rotation_test){
  PolarSite ps(1,"ps2");
  
  Eigen::Matrix3d R; //Rotation around z axes
  R << 0, -1, 0,
      1,  0,  0,
      0,  0,  1 ;
  
  Eigen::Vector3d shift;
  shift<<0,0,1;
  ps.Translate(shift);
  Eigen::Vector3d pos=ps.getPos(); // pos=(0,0,1)



  Eigen::VectorXd multipoles=Eigen::VectorXd::Zero(9);
  multipoles<<1,1,0,0,0,1,0,0,0; //q=1, mu_x=1 and Q_21c=1 the rest is 0
  
  ps.setMultipoles(multipoles);
  
  
Eigen::VectorXd rotmultipoles=Eigen::VectorXd::Zero(9);
rotmultipoles<<1,0,1,0,0,0,1,0,0; //q=1, mu_y=1 and Q_21s=1 is 0
  

ps.Rotate(R);
  
bool equalpos=pos.isApprox(ps.getPos(),1e-5);

bool equalmultipoles=rotmultipoles.isApprox(ps.getMultipoles(),1e-5);


  BOOST_CHECK_EQUAL(equalpos,true); //pos not affected by rotation
  BOOST_CHECK_EQUAL(equalmultipoles,true); 
}

BOOST_AUTO_TEST_CASE(interaction_test) {
  PolarSite ps1(1,"ps1");
  PolarSite ps2(2,"ps2");
    
  Eigen::Vector3d shift;
  shift<<1,0,0;
  ps2.Translate(shift);
  
  
  Eigen::VectorXd mp1 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd mp2 = Eigen::VectorXd::Zero(1);
  
  mp1<<1;
  mp2<<-1;
  ps1.setMultipoles(mp1);
  ps2.setMultipoles(mp2);
  
 
  double Energyref=-1;
  double Energy= ps1.InteractionAB(ps2);
   BOOST_CHECK_EQUAL(Energy==Energyref,true); 
  
}



BOOST_AUTO_TEST_SUITE_END()
