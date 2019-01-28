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

#ifndef _VOTCA_XTP_BSE_OP_H
#define _VOTCA_XTP_BSE_OP_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/matrixfreeoperator.h>

namespace votca {
namespace xtp {

class BSE_OPERATOR : public MatrixFreeOperator {

 public:
 
  BSE_OPERATOR(Orbitals& orbitals,ctp::Logger &log,TCMatrix_gwbse& Mmn,const Eigen::MatrixXd& Hqp):
        _log(log),
        _orbitals(orbitals),
        //_eh_s(orbitals.eh_s()),
        //_eh_t(orbitals.eh_t()),
        //_bse_singlet_energies(orbitals.BSESingletEnergies()),
        //_bse_singlet_coefficients(orbitals.BSESingletCoefficients()),
        //_bse_singlet_coefficients_AR(orbitals.BSESingletCoefficientsAR()),
        //_bse_triplet_energies(orbitals.BSETripletEnergies()),
        //_bse_triplet_coefficients(orbitals.BSETripletCoefficients()),
        _Mmn(Mmn),_Hqp(Hqp){};

    struct options {
        bool useTDA=true;
        int homo;
        int rpamin;
        int rpamax;
        int qpmin;
        int vmin;
        int cmax;
        int nmax; //number of eigenvectors to calculate
        bool davidson=0; // use davidson to diagonalize the matrix
        bool jocc=0; // jacobi orthogonal correction instead of DPR
        int jocc_linsolve=1; //method to solve the linea system in jacobi davidson
        double min_print_weight=0.5;  //minimium contribution for state to print it
        };
  
 

   void configure(const options& opt){
    _opt=opt;
    _bse_vmax = _opt.homo;
    _bse_cmin = _opt.homo+1;
    _bse_vtotal = _bse_vmax - _opt.vmin + 1;
    _bse_ctotal =_opt.cmax - _bse_cmin + 1;
    _bse_size = _bse_vtotal * _bse_ctotal;
    SetupDirectInteractionOperator();
  }
  
    
 private:
    options _opt;

    struct Interaction {
      Eigen::VectorXd exchange_contrib;
      Eigen::VectorXd direct_contrib;
      Eigen::VectorXd qp_contrib;
    };
      
  ctp::Logger &_log;
  int  _bse_vmax;
  int  _bse_cmin;
  int  _bse_size;
  int  _bse_vtotal;
  int  _bse_ctotal;
  
  Orbitals& _orbitals;
    // BSE variables and functions
  //MatrixXfd& _eh_s;  // only for storage in orbitals object
  //MatrixXfd& _eh_t;  // only for storage in orbitals object
  
   // references are stored in orbitals object
  // VectorXfd& _bse_singlet_energies;  
  // MatrixXfd& _bse_singlet_coefficients;                                                 
  // MatrixXfd& _bse_singlet_coefficients_AR;  
  // VectorXfd& _bse_triplet_energies;  
  // MatrixXfd& _bse_triplet_coefficients; 
  
  TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _Hqp;

  VectorXfd _epsilon_0_inv;

  void SetupDirectInteractionOperator();

 protected :

  Eigen::VectorXd Hqp_col(int index) const;

  template <int factor>
  Eigen::VectorXd Hx_col(int index) const;

  Eigen::VectorXd Hd_col(int index) const;

  template <int factor>
  Eigen::VectorXd Hd2_col(int index) const;


 
 
};
}
}

#endif /* _VOTCA_XTP_BSE_OP_H */
