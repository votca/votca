/* 
 *           Copyright 2009-2018 The VOTCA Development Team
 *                      (http://www.votca.org)
 *
 *     Licensed under the Apache License,Version 2.0 (the "License")
 *
 *You may not use this file except in compliance with the License.
 *You may obtain a copy of the License at
 *
 *             http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing,software
 *distributed under the License is distributed on an "AS IS" BASIS,
 *WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,either express or implied.
 *See the License for the specific language governing permissions and
 *limitations under the License.
 *
 */


#include <votca/xtp/staticsite.h>
#include <boost/format.hpp>
#include <fstream>
#include <string>
#include <votca/tools/constants.h>

#include "votca/xtp/interaction.h"


using namespace std;

namespace votca {
  namespace xtp {


    template< class S, class T>
    Eigen::MatrixXd Interactor::FillInteraction(const S& Site1, T& Site2) {

      const Eigen::Vector3d& posB = Site2.getPos();
      // const Eigen::Vector3d& posA = this->getPos();
      const Eigen::Vector3d& posA = Site1.getPos();
      
      int rankA = Site1.getRank();
      int rankB = Site2.getRank();
      int size1 = Site1.getPermMultipole().size();
      if (Site1.Identify()=="polarsite" && rankA < 1) {
        rankA = 1;
        size1 = 4;
      }
      int size2 = Site2.getPermMultipole().size();
      if (Site2.Identify()=="polarsite" && rankB < 1) {
        rankB = 1;
        size2 = 4;
      }

      Eigen::MatrixXd interaction = Eigen::MatrixXd::Zero(size1, size2);
      const Eigen::Vector3d r_AB = posB - posA; //Vector of the distance between polar sites
      const double R = r_AB.norm(); //Norm of distance vector
      const Eigen::Vector3d pos_a = r_AB / R; //unit vector on the sites reciprocal direction; This points toward A
      const Eigen::Vector3d pos_b = -pos_a; //unit vector on the sites reciprocal direction; This points toward B   

      const double fac0 = 1 / R;

      //Charge-Charge Interaction
      interaction(0, 0) = fac0; //T_00,00
      const double sqr3 = std::sqrt(3);
      const Eigen::Matrix3d AxA = pos_a * pos_a.transpose();
      const Eigen::Matrix3d&BxB = AxA;

      const double fac1 = std::pow(fac0, 2);
      const double fac2 = std::pow(fac0, 3);

      if (rankA > 0) {
        //Dipole-Charge Interaction
        interaction.block<3, 1>(1, 0) = fac1*pos_a; //T_1alpha,00 (alpha=x,y,z)
      }
      if (rankA > 1) {
        //Quadrupole-Charge interaction
        interaction(4, 0) = fac2 * 0.5 * (3 * AxA(2, 2) - 1); //T20,00
        interaction(5, 0) = fac2 * sqr3 * AxA(0, 2); //21c,00
        interaction(6, 0) = fac2 * sqr3 * AxA(1, 2); //T21s,000
        interaction(7, 0) = fac2 * 0.5 * sqr3 * (AxA(0, 0) - AxA(1, 1)); //T22c,00
        interaction(8, 0) = fac2 * sqr3 * AxA(0, 1); //T22s,00
      }

      if (rankB > 0) {
        //Charge-Dipole Interaction
        interaction.block<1, 3>(0, 1) = fac1*pos_b; //T_00,1alpha (alpha=x,y,z)
      }
      if (rankB > 1) {
        //Charge-Quadrupole Interaction
        interaction(0, 4) = fac2 * 0.5 * (3 * BxB(2, 2) - 1); //T00,20
        interaction(0, 5) = fac2 * sqr3 * BxB(0, 2); //T00,21c
        interaction(0, 6) = fac2 * sqr3 * BxB(1, 2); //T00,21s
        interaction(0, 7) = fac2 * 0.5 * sqr3 * (BxB(0, 0) - BxB(1, 1)); //T00,22c
        interaction(0, 8) = fac2 * sqr3 * BxB(0, 1); //T00,22s
      }

      const Eigen::Matrix3d c = Eigen::Matrix3d::Identity();
      if (rankA > 0 && rankB > 0) {
        const Eigen::Matrix3d AxB = pos_a * pos_b.transpose();
        //Dipole-Dipole Interaction 
        interaction.block<3, 3>(1, 1) = fac2 * (3 * AxB + c); //T_1alpha,1beta (alpha,beta=x,y,z) 
        const double fac3 = std::pow(fac0, 4);
        if (rankA > 1) {
          //Quadrupole-Dipole Interaction
          interaction.block<1, 3>(4, 1) = 0.5 * fac3 * (15 * AxA(2, 2) * pos_b.transpose() + 6 * pos_a(2) * c.row(2) - 3 * pos_b.transpose()); //T20-1beta (beta=x,y,z)
          interaction.block<1, 3>(5, 1) = fac3 * sqr3 * (pos_a(0) * c.row(2) + c.row(0) * pos_a(2) + 5 * AxA(0, 2) * pos_b.transpose()); //T21c-1beta (beta=x,y,z)
          interaction.block<1, 3>(6, 1) = fac3 * sqr3 * (pos_a(1) * c.row(2) + c.row(1) * pos_a(2) + 5 * AxA(1, 2) * pos_a.transpose()); //T21s-1beta (beta=x,y,z)
          interaction.block<1, 3>(7, 1) = fac3 * 0.5 * sqr3 * (5 * (AxA(0, 0) - AxA(1, 1)) * pos_b.transpose() + 2 * pos_a(0) * c.row(0) - 2 * pos_a(1) * c.row(1)); //T22c-1beta (beta=x,y,z)
          interaction.block<1, 3>(8, 1) = fac3 * sqr3 * (5 * AxA(0, 1) * pos_b.transpose() + pos_a(0) * c.row(1) + pos_a(1) * c.row(0)); //T22s-1beta (beta=x,y,z)    
        }
        if (rankB > 1) {
          //Dipole-Quadrupole Interaction
          interaction.block<3, 1>(1, 4) = 0.5 * fac3 * (15 * BxB(2, 2) * pos_a + 6 * pos_b(2) * c.col(2) - 3 * pos_a); //T1beta-20 (beta=x,y,z)
          interaction.block<3, 1>(1, 5) = fac3 * sqr3 * (pos_b(0) * c.col(2) + c.col(0) * pos_b(2) + 5 * BxB(0, 2) * pos_a); //T1beta-21c (beta=x,y,z)
          interaction.block<3, 1>(1, 6) = fac3 * sqr3 * (pos_b(1) * c.col(2) + c.col(1) * pos_b(2) + 5 * BxB(1, 2) * pos_b); //T1beta-21s (beta=x,y,z)
          interaction.block<3, 1>(1, 7) = 0.5 * fac3 * sqr3 * (5 * (BxB(0, 0) - BxB(1, 1)) * pos_a + 2 * pos_b(0) * c.col(0) - 2 * pos_b(1) * c.col(1)); //T1beta-22c (beta=x,y,z)
          interaction.block<3, 1>(1, 8) = fac3 * sqr3 * (5 * BxB(0, 1) * pos_a + pos_b(0) * c.col(1) + pos_b(1) * c.col(0)); //T1beta-22s (beta=x,y,z)    
        }

        if (rankA > 1 && rankB > 1) {
          const double fac4 = std::pow(fac0, 5);
          //Quadrupole-Quadrupole Interaction
          interaction(4, 4) = fac4 * (3. / 4.)*(35 * AxA(2, 2) * BxB(2, 2) - 5 * AxA(2, 2) - 5 * BxB(2, 2)
                  + 20 * AxB(2, 2) * c(2, 2) + 2 * c(2, 2) * c(2, 2) + 1); //T20,20
          interaction(4, 5) = 0.5 * fac4 * sqr3 * (35 * AxA(2, 2) * BxB(0, 2) - 5 * BxB(0, 2) + 10 * AxB(2, 0) * c(2, 2)
                  + 10 * AxB(2, 2) * c(2, 1) + 2 * c(2, 0) * c(2, 2)); //T20,21c
          interaction(4, 6) = 0.5 * fac4 * sqr3 * (35 * AxA(2, 2) * BxB(1, 2) - 5 * BxB(1, 2) + 10 * AxB(2, 1) * c(2, 2)
                  + 10 * AxB(2, 2) * c(2, 1) + 2 * c(2, 1) * c(2, 2)); //T20,21s
          interaction(4, 7) = 0.25 * fac4 * sqr3 * (35 * AxB(2, 0) - 35 * AxB(2, 1) - 5 * BxB(0, 0)
                  + 5 * BxB(1, 1) + 20 * AxB(2, 0) * c(2, 0) - 20 * AxB(2, 1) * c(2, 1)
                  + 2 * c(2, 0) * c(2, 0) - 2 * c(2, 1) * c(2, 1)); //T20,22c
          interaction(4, 8) = 0.5 * fac4 * sqr3 * (35 * AxA(2, 2) * BxB(0, 1) - 5 * BxB(0, 1)
                  + 10 * AxB(2, 0) * c(2, 1) + 10 * AxB(2, 1) * c(2, 0) + 2 * c(2, 0) * c(2, 1)); //T20,22s
          interaction(5, 5) = fac4 * (35 * AxA(0, 2) * BxB(0, 2) + 5 * AxB(0, 0) * c(2, 2) + 5 * AxB(0, 2) * c(2, 0)
                  + 5 * AxB(2, 0) * c(0, 2) + 5 * AxB(2, 2) * c(0, 0) + c(0, 0) * c(2, 2) + c(0, 2) * c(2, 0)); //T21c,21c
          interaction(5, 6) = fac4 * (35 * AxA(0, 2) * BxB(1, 2) + 5 * AxB(0, 1) * c(2, 2) + 5 * AxB(0, 2) * c(2, 1)
                  + 5 * AxB(2, 1) * c(0, 2) + 5 * AxB(2, 2) * c(0, 1) + c(0, 1) * c(2, 2) + c(0, 2) * c(2, 1)); //T21c,21s
          interaction(5, 7) = 0.5 * fac4 * (35 * AxA(0, 2) * BxB(0, 0) - 35 * AxA(0, 2) * BxB(1, 1) + 10 * AxB(0, 0) * c(2, 0)
                  - 10 * AxB(0, 1) * c(2, 1) + 10 * AxB(0, 0) * c(0, 0) - 10 * AxB(2, 1) * c(0, 1) + 2 * c(0, 0) * c(2, 0) - 2 * c(0, 1) * c(2, 1)); //T21c,22c
          interaction(5, 8) = fac4 * (35 * AxA(0, 2) * BxB(0, 1) + 5 * AxB(0, 0) * c(2, 1) + 5 * AxB(0, 1) * c(2, 0)
                  + 5 * AxB(2, 0) * c(0, 1) + 5 * AxB(2, 1) * c(0, 0) + c(0, 0) * c(2, 1) + c(0, 1) * c(2, 0)); //T21c,22s
          interaction(6, 6) = fac4 * (35 * AxA(1, 2) * BxB(1, 2) + 5 * AxB(1, 1) * c(2, 2) + 5 * AxB(1, 2) * c(2, 1)
                  + 5 * AxB(2, 1) * c(1, 2) + 5 * AxB(2, 2) * c(1, 1) + c(1, 1) * c(2, 2) + c(1, 2) * c(2, 1)); //T21s,21s
          interaction(6, 7) = 0.5 * fac4 * (35 * AxA(1, 2) * BxB(0, 0) - 35 * AxA(1, 2) * BxB(1, 1) + 10 * AxB(1, 2) * c(2, 0)
                  - 10 * AxB(1, 1) * c(2, 1) + 10 * AxB(2, 0) * c(1, 0) - 10 * AxB(2, 1) * c(1, 1) + 2 * c(1, 0) * c(2, 0) - 2 * c(1, 1) * c(2, 1)); //T21s,22c
          interaction(6, 8) = fac4 * (35 * AxA(1, 2) * BxB(0, 1) + 5 * AxB(1, 0) * c(2, 1) + 5 * AxB(1, 1) * c(2, 1)
                  + 5 * AxB(2, 0) * c(1, 1) + 5 * AxB(2, 1) * c(1, 2) + c(1, 0) * c(2, 1) + c(1, 1) * c(2, 0)); //T21s,22s
          interaction(7, 7) = 0.25 * fac4 * (35 * AxA(0, 0) * BxB(0, 0) - 35 * AxA(0, 0) * BxB(1, 1)
                  - 35 * AxA(1, 1) * BxB(0, 0) + 35 * AxA(1, 1) * BxB(1, 1) + 20 * AxB(0, 0) * c(0, 0)
                  - 20 * AxB(0, 1) * c(0, 1) - 20 * AxB(1, 0) * c(1, 0) + 20 * AxB(0, 0) * c(1, 1)
                  + 2 * c(0, 0) * c(0, 0) - 2 * c(0, 1) * c(0, 1) - 2 * c(1, 0) * c(1, 0) + 2 * c(1, 1) * c(1, 1)); //T22c,22c
          interaction(7, 8) = 0.5 * fac4 * (35 * BxB(0, 1) * AxA(0, 0) - 35 * BxB(1, 2) * AxA(1, 1)
                  + 10 * AxB(0, 0) * c(0, 1) + 10 * AxB(0, 1) * c(0, 0) - 10 * AxB(1, 0) * c(1, 1)
                  - 10 * AxB(1, 1) * c(1, 2) + 2 * c(0, 0) * c(0, 1) - 2 * c(1, 0) * c(1, 1)); //T22c,22s
          interaction(8, 8) = 0.5 * fac4 * (35 * AxA(0, 1) * BxB(0, 1) + 5 * AxB(0, 0) * c(1, 1)
                  + 5 * AxB(0, 1) * c(1, 0) + 5 * AxB(1, 0) * c(0, 1) + 5 * AxB(1, 1) * c(0, 0) + c(0, 0) * c(1, 1) + c(0, 1) * c(1, 0)); //T22s,22s
          interaction(5, 4) = 0.5 * fac4 * sqr3 * (35 * BxB(2, 2) * AxA(0, 2) - 5 * AxA(0, 2)
                  + 10 * AxB(0, 2) * c(2, 2) + 10 * AxB(2, 2) * c(1, 2) + 2 * c(0, 2) * c(2, 2)); //T21c,20
          interaction(6, 4) = 0.5 * fac4 * sqr3 * (35 * BxB(2, 2) * AxA(1, 2) - 5 * AxA(1, 2)
                  + 10 * AxB(1, 2) * c(2, 2) + 10 * AxB(2, 2) * c(1, 2) + 2 * c(1, 2) * c(2, 2)); //T21s,20
          interaction(6, 5) = fac4 * (35 * BxB(0, 2) * AxA(1, 2) + 5 * AxB(1, 0) * c(2, 2)
                  + 5 * AxB(2, 0) * c(1, 2) + 5 * AxB(2, 1) * c(2, 0) + 5 * AxB(2, 2) * c(1, 0) + c(1, 0) * c(2, 2) + c(2, 0) * c(1, 2)); //T21s,21c
          interaction(7, 4) = 0.25 * fac4 * sqr3 * (35 * AxB(0, 2) - 35 * BxB(2, 2) * AxA(1, 1)
                  - 5 * AxA(0, 0) + 5 * AxA(1, 1) + 20 * AxB(0, 2) * c(0, 2) - 20 * AxB(1, 2) * c(1, 2)
                  + 2 * c(0, 2) * c(0, 2) - 2 * c(1, 2) * c(1, 2)); //T22c,20
          interaction(7, 5) = 0.5 * fac4 * (35 * BxB(0, 2) * AxA(0, 0) - 35 * BxB(0, 2) * AxA(1, 1)
                  + 10 * AxB(0, 0) * c(0, 2) - 10 * AxB(1, 0) * c(1, 2) + 10 * AxB(0, 0) * c(0, 0)
                  - 10 * AxB(1, 2) * c(1, 0) + 2 * c(0, 0) * c(0, 2) - 2 * c(1, 0) * c(1, 2)); //T22c,21c
          interaction(7, 6) = 0.5 * fac4 * (35 * BxB(1, 2) * AxA(0, 0) - 35 * BxB(1, 2) * AxA(1, 1)
                  + 10 * AxB(2, 1) * c(0, 2) - 10 * AxB(1, 1) * c(1, 2) + 10 * AxB(0, 2) * c(0, 1)
                  - 10 * AxB(1, 2) * c(1, 1) + 2 * c(0, 1) * c(0, 2) - 2 * c(1, 1) * c(1, 2)); //T22c,21s
          interaction(8, 4) = 0.5 * fac4 * sqr3 * (35 * BxB(2, 2) * AxA(0, 1) - 5 * AxA(0, 1)
                  + 10 * AxB(0, 2) * c(1, 2) + 10 * AxB(1, 2) * c(0, 2) + 2 * c(0, 2) * c(1, 2)); //T22s,20
          interaction(8, 5) = fac4 * (35 * BxB(0, 2) * AxA(1, 0) + 5 * AxB(0, 0) * c(1, 2)
                  + 5 * AxB(1, 0) * c(0, 2) + 5 * AxB(0, 2) * c(1, 0) + 5 * AxB(1, 2) * c(0, 0)
                  + c(0, 0) * c(1, 2) + c(1, 0) * c(0, 2)); //T22s,21c
          interaction(8, 6) = fac4 * (35 * BxB(1, 2) * AxA(0, 1) + 5 * AxB(0, 1) * c(1, 2)
                  + 5 * AxB(1, 1) * c(1, 2) + 5 * AxB(0, 2) * c(1, 1) + 5 * AxB(1, 2) * c(2, 1)
                  + c(0, 1) * c(1, 2) + c(1, 1) * c(0, 2)); //T22s,21s
          interaction(8, 7) = 0.5 * fac4 * (35 * AxA(1, 0) * BxB(0, 0) - 35 * AxA(1, 2) * BxB(1, 1)
                  + 10 * AxB(0, 0) * c(1, 0) + 10 * AxB(1, 0) * c(0, 0) - 10 * AxB(0, 1) * c(1, 1)
                  - 10 * AxB(1, 1) * c(2, 1) + 2 * c(0, 0) * c(1, 0) - 2 * c(0, 1) * c(1, 1)); //T22s,22c          
        }

      }
      return interaction; //in units of 4piepsilon0  
    }

    Eigen::MatrixXd Interaction::FillTholeInteraction(const S& Site1, T& Site2, double a) {

      const Eigen::Vector3d& posB = Site2.getPos();
      const Eigen::Vector3d& posA = Site1.getPos();
      int rankA = Site1.getRank();
      int rankB = Site2.getRank();
      int size1 = Site1.getPermMultipole().size();
      if (rankA < 1) {
        rankA = 1;
        size1 = 4;
      }
      int size2 = Site2.getPermMultipole().size();
      if (rankB < 1) {
        rankB = 1;
        size2 = 4;
      }

      Eigen::MatrixXd interaction = Eigen::MatrixXd::Zero(size1, size2);
      const Eigen::Vector3d r_AB = posB - posA; //Vector of the distance between polar sites
      const double R = r_AB.norm(); //Norm of distance vector
      const Eigen::Vector3d pos_a = r_AB / R; //unit vector on the sites reciprocal direction; This points toward A
      const Eigen::Vector3d pos_b = -pos_a; //unit vector on the sites reciprocal direction; This points toward B   

      const double sqr3 = std::sqrt(3);
      const Eigen::Matrix3d AxA = pos_a * pos_a.transpose();
      const Eigen::Matrix3d& BxB = AxA;
      const double fac0 = 1 / R;
      const double fac1 = std::pow(fac0, 2);
      const double fac2 = std::pow(fac0, 3);
      const double au3 = a / (fac2 * std::sqrt(Site1._eigendamp * Site2._eigendamp));
      double lambda3 = 1.0;
      double lambda5 = 1.0;
      double lambda7 = 1.0;
      if (au3 < 40) {
        const double exp_ua = std::exp(-au3);
        lambda3 = 1 - exp_ua;
        lambda5 = 1 - (1 + au3) * exp_ua;
        lambda7 = 1 - (1 + au3 + 0.6 * au3 * au3) * exp_ua;
      }
      if (rankA > 0) {
        //Dipole-Charge Interaction
        interaction.block<3, 1>(1, 0) = lambda3 * fac1*pos_a; //T_1alpha,00 (alpha=x,y,z)
      }
      if (rankB > 0) {
        //Charge-Dipole Interaction
        interaction.block<1, 3>(0, 1) = lambda3 * fac1*pos_b; //T_00,1alpha (alpha=x,y,z)
      }

      const Eigen::Matrix3d c = Eigen::Matrix3d::Identity();
      if (rankA > 0 && rankB > 0) {
        const Eigen::Matrix3d AxB = pos_a * pos_b.transpose();
        //Dipole-Dipole Interaction 
        interaction.block<3, 3>(1, 1) = fac2 * (lambda5 * 3 * AxB + lambda3 * c); //T_1alpha,1beta (alpha,beta=x,y,z) 
        const double fac3 = std::pow(fac0, 4);
        if (rankA > 1) {
          //Quadrupole-Dipole Interaction
          interaction.block<1, 3>(4, 1) = 0.5 * fac3 * (lambda7 * 15 * AxA(2, 2) * pos_b.transpose() + lambda5 * (6 * pos_a(2) * c.row(2) - 3 * pos_b.transpose())); //T20-1beta (beta=x,y,z)
          interaction.block<1, 3>(5, 1) = fac3 * sqr3 * (lambda5 * (pos_a(0) * c.row(2) + c.row(0) * pos_a(2)) + lambda7 * 5 * AxA(0, 2) * pos_b.transpose()); //T21c-1beta (beta=x,y,z)
          interaction.block<1, 3>(6, 1) = fac3 * sqr3 * (lambda5 * (pos_a(1) * c.row(2) + c.row(1) * pos_a(2)) + lambda7 * 5 * AxA(1, 2) * pos_a.transpose()); //T21s-1beta (beta=x,y,z)
          interaction.block<1, 3>(7, 1) = fac3 * 0.5 * sqr3 * (lambda7 * 5 * (AxA(0, 0) - AxA(1, 1)) * pos_b.transpose() + lambda5 * (2 * pos_a(0) * c.row(0) - 2 * pos_a(1) * c.row(1))); //T22c-1beta (beta=x,y,z)
          interaction.block<1, 3>(8, 1) = fac3 * sqr3 * (lambda7 * 5 * AxA(0, 1) * pos_b.transpose() + lambda5 * (pos_a(0) * c.row(1) + pos_a(1) * c.row(0))); //T22s-1beta (beta=x,y,z)    
        }
        if (rankB > 1) {
          //Dipole-Quadrupole Interaction
          interaction.block<3, 1>(1, 4) = 0.5 * fac3 * (lambda7 * 15 * BxB(2, 2) * pos_a + lambda5 * (6 * pos_b(2) * c.col(2) - 3 * pos_a)); //T1beta-20 (beta=x,y,z)
          interaction.block<3, 1>(1, 5) = fac3 * sqr3 * (lambda5 * (pos_b(0) * c.col(2) + c.col(0) * pos_b(2)) + lambda7 * 5 * BxB(0, 2) * pos_a); //T1beta-21c (beta=x,y,z)
          interaction.block<3, 1>(1, 6) = fac3 * sqr3 * (lambda5 * (pos_b(1) * c.col(2) + c.col(1) * pos_b(2)) + lambda7 * 5 * BxB(1, 2) * pos_b); //T1beta-21s (beta=x,y,z)
          interaction.block<3, 1>(1, 7) = 0.5 * fac3 * sqr3 * (lambda7 * 5 * (BxB(0, 0) - BxB(1, 1)) * pos_a + lambda5 * (2 * pos_b(0) * c.col(0) - 2 * pos_b(1) * c.col(1))); //T1beta-22c (beta=x,y,z)
          interaction.block<3, 1>(1, 8) = fac3 * sqr3 * (lambda7 * 5 * BxB(0, 1) * pos_a + lambda5 * (pos_b(0) * c.col(1) + pos_b(1) * c.col(0))); //T1beta-22s (beta=x,y,z)    
        }

      }
      return interaction; //in units of 4piepsilon0  
    }
    
    double Interactor::InteractStatic(StaticSite& Site1, StaticSite& Site2){
      const Eigen::MatrixXd Tab = FillInteraction(Site1,Site2); // T^(ab)_tu
      const Eigen::VectorXd& multipolesA = Site1.getPermMultipole(); //Q^(a)_t
      const Eigen::VectorXd& multipolesB = Site2.getPermMultipole(); //Q^(b)_u  
      
      int sizeB=multipolesB.size();
      int sizeA=multipolesA.size();
      //Potentials
      this->PhiP += (Tab.row(0).head(sizeB) * multipolesB).sum();
      otherSite.PhiP += (multipolesA.transpose() * Tab.col(0).head(sizeA)).sum();
      //Fields
      if (this->isPolarisable()) {
        this->_localpermanetField += (Tab.block(1, 0, 3, sizeB) * multipolesB).rowwise().sum();
      }
      if (Site2.isPolarisable()) {
        otherSite._localpermanetField += (multipolesA.transpose() * Tab.block(0, 1, sizeA, 3)).colwise().sum().transpose();
      }
      //Energy
      const double EnergyAB = multipolesA.transpose() * Tab.topLeftCorner(sizeA,sizeB)*multipolesB; //Interaction Energy between sites A and B

      return EnergyAB;
    }
    
    double Interactor::InteractInduction(PolarSite& Site1, PolarSite& Site2, double a = 0.39){
      const Eigen::VectorXd& multipolesA = Site1.getPermMultipole(); //Q^(a)_t
    
      const Eigen::VectorXd& multipolesB = Site2.getPermMultipole(); //Q^(b)_u  
      double EnergyAB = 0;
      const Eigen::MatrixXd tTab = FillTholeInteraction(Site2, a); // \tilde{T}^(ab)_tu
      //Calculate Potential due to induced Dipoles
      if (Site2.isPolarisable()) {
        this->PhiU += (tTab.row(0).segment<3>(1)) * Site2._inducedDipole;
        EnergyAB += multipolesA.transpose() * tTab.block(0, 1, multipolesA.size(), 3) * otherSite._inducedDipole;
      }
      if (this->isPolarisable()) {
        otherSite.PhiP += (this->_inducedDipole.transpose() * tTab.col(0).segment<3>(1)).sum();
        EnergyAB += this->_inducedDipole.transpose() * tTab.block(1, 0, 3, multipolesB.size()) * multipolesB;
        //Calculate Interaction induced Dipole static Poles
      }
      if (this->isPolarisable() && otherSite.isPolarisable()) {
        //Calculate Field due to induced Dipoles
        Eigen::Matrix3d Tdd = tTab.block<3, 3>(1, 1);
        this->_localinducedField += (Tdd * otherSite._inducedDipole).rowwise().sum();
        otherSite._localinducedField += (this->_inducedDipole.transpose() * Tdd).colwise().sum().transpose();
        //Calculate Interaction induced Dipole induced Dipole
        EnergyAB += this->_inducedDipole.transpose() * Tdd * otherSite._inducedDipole;

      }
      return EnergyAB;
    }


  }
}
