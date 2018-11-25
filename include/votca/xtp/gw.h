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

#ifndef _VOTCA_XTP_GW_H
#define _VOTCA_XTP_GW_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

class GW {     
 public:
    GW(ctp::Logger &log, TCMatrix_gwbse& Mmn):_log(log),_Mmn(Mmn){};

    void configure(int homo, int qpmin,int qpmax,int g_sc_max_iterations,
                double g_sc_limit){
      _homo=homo;
      _qpmin=qpmin;
      _qpmax=qpmax;
      _qptotal=_qpmax - _qpmin + 1;
      _g_sc_limit=g_sc_limit;
      _g_sc_max_iterations=g_sc_max_iterations;
      if(_g_sc_max_iterations<1) {_g_sc_max_iterations=1;}
  }

    void setSelfConsistencyOptions



 
 
 private:

    ctp::Logger &_log;
    TCMatrix_gwbse& _Mmn;
     
 void PrintQP_Energies(const Eigen::VectorXd& gwa_energies, const Eigen::VectorXd& qp_diag_energies)const;
 void PrintGWA_Energies(const Eigen::MatrixXd& vxc,const Sigma& sigma, const Eigen::VectorXd& dft_energies)const;
      

};
}
}

#endif /* _VOTCA_XTP_BSE_H */
