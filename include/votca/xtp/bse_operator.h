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

#ifndef _VOTCA_XTP_BSE_OPERATOR_H
#define _VOTCA_XTP_BSE_OPERATOR_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/matrixfreeoperator.h>
//#include <votca/xtp/bse_engine.h>

namespace votca {
namespace xtp {

class BSE_OPERATOR : public MatrixFreeOperator {

 public:
 
  BSE_OPERATOR(Orbitals& orbitals,ctp::Logger &log,TCMatrix_gwbse& Mmn,const Eigen::MatrixXd& Hqp):
        _log(log),
        _orbitals(orbitals),
        _Mmn(Mmn),_Hqp(Hqp){};

    struct options {
        int homo; //
        int rpamin; //
        int rpamax; //
        int qpmin;  //
        int vmin; // 
        int cmax; // 
        };
  
  options _opt;

  int  _bse_vmax;
  int  _bse_cmin;
  int  _bse_size;
  int  _bse_vtotal;
  int  _bse_ctotal;    

  void SetupDirectInteractionOperator();

 private:

    struct Interaction {
      Eigen::VectorXd exchange_contrib;
      Eigen::VectorXd direct_contrib;
      Eigen::VectorXd qp_contrib;
    };
      
  ctp::Logger &_log;

  Orbitals& _orbitals;
  
  TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _Hqp;

  VectorXfd _epsilon_0_inv;



  protected: 

  Eigen::VectorXd Hqp_col(int index) const;

  //template <int factor>
  Eigen::VectorXd Hx_col(int index) const;

  Eigen::VectorXd Hd_col(int index) const;

  //template <int factor>
  Eigen::VectorXd Hd2_col(int index) const;


 
 
};
}
}

#endif /* _VOTCA_XTP_BSE_OP_H */
