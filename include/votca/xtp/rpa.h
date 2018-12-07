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

#ifndef _VOTCA_XTP_RPA_H
#define _VOTCA_XTP_RPA_H
#include <votca/xtp/eigen.h>
#include <vector>


namespace votca
{
namespace xtp
{
class TCMatrix_gwbse;

class RPA
{
public:

    void configure(int homo, int rpamin, int rpamax){
        _homo = homo;
        _rpamin = rpamin;
        _rpamax = rpamax;
    }

    Eigen::MatrixXd calculate_epsilon_i(double frequency, const Eigen::VectorXd& qp_energies, const TCMatrix_gwbse& Mmn)const{
        return calculate_epsilon<true>(frequency, qp_energies, Mmn);
    }

    Eigen::MatrixXd calculate_epsilon_r(double frequency, const Eigen::VectorXd& qp_energies, const TCMatrix_gwbse& Mmn)const{
        return calculate_epsilon<false>(frequency, qp_energies, Mmn);
    }

private:

    int _homo; // HOMO index
    int _rpamin;
    int _rpamax;

    template< bool imag>
    Eigen::MatrixXd calculate_epsilon(double frequency, const Eigen::VectorXd& qp_energies, const TCMatrix_gwbse& Mmn)const;

};
}
}

#endif /* _VOTCA_RPA_RPA_H */
