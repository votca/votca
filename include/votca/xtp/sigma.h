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
#include <votca/ctp/logger.h>

namespace votca {
namespace xtp {
    class PPM;
    class TCMatrix_gwbse;

class Sigma {
 public:
  Sigma(const TCMatrix_gwbse& Mmn):_Mmn(Mmn){};
  
  void configure(int homo, int qpmin,int qpmax){
      _homo=homo;
      _qpmin=qpmin;
      _qpmax=qpmax;
      _qptotal=_qpmax - _qpmin + 1;
  }
  
Eigen::MatrixXd CalcExchange()const;

Eigen::VectorXd CalcCorrelationDiag(const PPM & ppm,const Eigen::VectorXd& energies )const;

Eigen::MatrixXd CalcCorrelationOffDiag(const PPM & ppm,const Eigen::VectorXd& energies)const;
 

 private:

  const TCMatrix_gwbse& _Mmn;

  inline double SumSymmetric(real_gwbse Mmn1xMmn2, double qpmin1, double qpmin2, const double gwa_energy);
  inline double Stabilize(double denom);
  inline void Stabilize(Eigen::ArrayXd& denom);
  int _homo;   // HOMO index
  int _qpmin;
  int _qpmax;
  int _qptotal;

 

  
};
}
}

#endif /* _VOTCA_XTP_SIGMA_H */
