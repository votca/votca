/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_XTP_SIGMA_H
#define _VOTCA_XTP_SIGMA_H
#include <votca/xtp/eigen.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {


class Sigma {
 public:
  Sigma(xtp::Logger *log){
_log = log;
_gwa_energies.resize(0);
  }
  
  void configure(unsigned homo, unsigned qpmin,unsigned qpmax,unsigned g_sc_max_iterations,
                double g_sc_limit){
      _homo=homo;
      _qpmin=qpmin;
      _qpmax=qpmax;
      _qptotal=_qpmax - _qpmin + 1;
      _g_sc_limit=g_sc_limit;
      _g_sc_max_iterations=g_sc_max_iterations;
      if(_g_sc_max_iterations<1) {_g_sc_max_iterations=1;}
  }
  
  void setDFTdata(double ScaHFX, const Eigen::MatrixXd* vxc, const Eigen::VectorXd* dftenergies){
      _ScaHFX=ScaHFX;
      _vxc=vxc;
      _dftenergies=dftenergies;
  }
  
void CalcdiagElements(const TCMatrix_gwbse& _Mmn,const PPM & ppm );
void CalcOffDiagElements(const TCMatrix_gwbse& _Mmn,const PPM & ppm );

Eigen::MatrixXd SetupFullQPHamiltonian();

const Eigen::VectorXd& getGWAEnergies()const{return _gwa_energies;}

void setGWAEnergies(const Eigen::VectorXd& gwa_energies){_gwa_energies=gwa_energies;}

double x(int i)const{return _sigma_x(i,i);}
double c(int i)const{return _sigma_c(i,i);}
  
void FreeMatrices(){
    _sigma_x.resize(0,0);
    _sigma_c.resize(0,0);
}

 private:
  void C_offdiag(const TCMatrix_gwbse& _Mmn, const PPM& ppm);
  void X_offdiag(const TCMatrix_gwbse& _Mmn);   
xtp::Logger *_log;
  unsigned _homo;   // HOMO index
  unsigned _qpmin;
  unsigned _qpmax;
  unsigned _qptotal;

  double _g_sc_limit;  // convergence criteria for g iteration [Hartree]]
  unsigned int _g_sc_max_iterations;
  
  double _ScaHFX;
  // Sigma related variables and functions
  Eigen::MatrixXd _sigma_x;  // exchange term
  Eigen::MatrixXd _sigma_c;  // correlation term

  
const Eigen::MatrixXd* _vxc;
const Eigen::VectorXd* _dftenergies;
  

  // QP variables and functions
  Eigen::VectorXd _gwa_energies;
  
};
}
}

#endif /* _VOTCA_XTP_SIGMA_H */
