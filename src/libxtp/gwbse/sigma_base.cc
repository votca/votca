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



#include <votca/xtp/sigma_base.h>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>

#include "votca/xtp/threecenter.h"


namespace votca {
  namespace xtp {

  Eigen::MatrixXd Sigma_base::CalcExchange()const{

    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(_qptotal,_qptotal);
    int gwsize = _Mmn.auxsize();
      #pragma omp parallel for schedule(dynamic)
      for (int gw_level1 = 0; gw_level1 < _qptotal; gw_level1++) {
        const MatrixXfd& Mmn1 = _Mmn[ gw_level1 + _qpmin ];
        for (int gw_level2 = gw_level1; gw_level2 < _qptotal; gw_level2++) {
          const MatrixXfd & Mmn2 = _Mmn[ gw_level2 + _qpmin ];
          double sigma_x =- (Mmn1.block(0,0,_homo+1,gwsize).cwiseProduct(Mmn2.block(0,0,_homo+1,gwsize))).sum();
          result(gw_level1, gw_level2) =  sigma_x;
          result(gw_level2, gw_level1) =  sigma_x;
        }
      }
        return result;
    }

       
  }
};
