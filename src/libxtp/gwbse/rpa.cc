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



#include <votca/xtp/rpa.h>
#include <votca/xtp/aomatrix.h>
#include "votca/xtp/threecenter.h"



namespace votca {
  namespace xtp {

 template< bool imag>
    Eigen::MatrixXd RPA::calculate_epsilon(double frequency)const{
        const int size = _Mmn.auxsize(); // size of gwbasis

        Eigen::MatrixXd result = Eigen::MatrixXd::Identity(size, size);
        const int lumo = _homo + 1;
        const int n_occ = lumo - _rpamin;
        const int n_unocc = _rpamax - _homo;
        const double freq2 = frequency*frequency;

#pragma omp parallel for
        for (int m_level = 0; m_level < n_occ; m_level++)        {
            const double qp_energy_m = _energies(m_level + _rpamin);
#if (GWBSE_DOUBLE)
            const Eigen::MatrixXd Mmn_RPA = _Mmn[ m_level].block(n_occ, 0, n_unocc, size);
#else
            const Eigen::MatrixXd Mmn_RPA = _Mmn[ m_level].block(n_occ, 0, n_unocc, size).cast<double>();
#endif
            const Eigen::ArrayXd deltaE=_energies.segment(lumo,n_unocc).array()-qp_energy_m;
            Eigen::VectorXd denom;
            if (imag){
                denom=4*deltaE/(deltaE.square()+freq2);
            }else{
                denom=2.0*((deltaE-frequency).inverse()+(deltaE+frequency).inverse());
            }
            Eigen::MatrixXd tempresult = Mmn_RPA.transpose() *denom.asDiagonal()* Mmn_RPA;

#pragma omp critical
            {
                result += tempresult;
            }
        }


        return result;
    }


 template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency)const;
 template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency)const;


  }}
