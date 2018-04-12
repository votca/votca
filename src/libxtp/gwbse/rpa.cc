/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

 void RPA::calculate_epsilon(const Eigen::VectorXd& qp_energies) {

            for (auto& matrix : _epsilon_r) {
                matrix = Eigen::MatrixXd::Identity(_Mmn_RPA.getAuxDimension(), _Mmn_RPA.getAuxDimension());
            }
            for (auto& matrix : _epsilon_i) {
                matrix = Eigen::MatrixXd::Identity(_Mmn_RPA.getAuxDimension(), _Mmn_RPA.getAuxDimension());
            }


            const int _size = _Mmn_RPA.getAuxDimension(); // size of gwbasis
            const int index_n = _Mmn_RPA.get_nmin();
            const int index_m = _Mmn_RPA.get_mmin();
#pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
                const double _qp_energy_m = qp_energies(_m_level + index_m);
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const Eigen::MatrixXd Mmn_RPA = _Mmn_RPA[ _m_level ].cast<double>();

#endif
                for (int i = 0; i < screen_freq_i.size(); ++i) {
                    // a temporary matrix, that will get filled in empty levels loop
                    double screen_freq2 = screen_freq_i(i) * screen_freq_i(i);
                    Eigen::VectorXd energydenom = Eigen::VectorXd::Zero(_Mmn_RPA.get_ntot());
                    for (int _n_level = 0; _n_level < energydenom.size(); _n_level++) {
                        const double _deltaE = qp_energies(_n_level + index_n) - _qp_energy_m;
                        energydenom(_n_level) = 4.0 * _deltaE / (_deltaE * _deltaE + screen_freq2); //hartree       
                    }
                    Eigen::MatrixXd tempresult = Mmn_RPA * energydenom.asDiagonal() * Mmn_RPA.transpose();

#pragma omp critical
                    {
                        _epsilon_i[i] += tempresult;
                    }
                }

                //real parts
                for (int i = 0; i < screen_freq_r.size(); ++i) {
                    Eigen::VectorXd energydenom = Eigen::VectorXd::Zero(_Mmn_RPA.get_ntot());
                    for (int _n_level = 0; _n_level < energydenom.size(); _n_level++) {
                        const double _deltaE = qp_energies(_n_level + index_n) - _qp_energy_m;
                        energydenom(_n_level) = 2.0 * (1.0 / (_deltaE - screen_freq_r(i)) + 1.0 / (_deltaE + screen_freq_r(i))); //hartree
                    }
                    Eigen::MatrixXd tempresult = Mmn_RPA * energydenom.asDiagonal() * Mmn_RPA.transpose();

#pragma omp critical
                    {
                        _epsilon_r[i] += tempresult;
                    }
                }

            } // occupied levels


            return;
        }

    void RPA::prepare_threecenters(const TCMatrix_gwbse& _Mmn_full) {
      //TODO maybe remove and instead not make a copy but use full Mmn with restricted indices
      _Mmn_RPA.Initialize(_Mmn_full.getAuxDimension(), _rpamin, _homo, _homo + 1,
              _rpamax);
      unsigned start = _Mmn_RPA.get_nmin() - _Mmn_full.get_nmin();
      // loop over m-levels in _Mmn_RPA
#pragma omp parallel for 
      for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
        // copy to _Mmn_RPA
        _Mmn_RPA[ _m_level ] = _Mmn_full[ _m_level ].block(0, start, _Mmn_full.getAuxDimension(), _Mmn_RPA.get_ntot());
      }// loop m-levels
      return;
    }





  }
};
