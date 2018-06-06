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



#include <votca/xtp/sigma.h>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>



namespace votca {
  namespace xtp {

    Eigen::MatrixXd Sigma::SetupFullQPHamiltonian() {

      // constructing full QP Hamiltonian
      Eigen::MatrixXd Hqp = _sigma_x + _sigma_c - (*_vxc);
      // diagonal elements are given by _qp_energies
      for (unsigned _m = 0; _m < Hqp.rows(); _m++) {
        Hqp(_m, _m) = _gwa_energies(_m + _qpmin);
      }
      return Hqp;
    }

    void Sigma::CalcdiagElements(const TCMatrix_gwbse& _Mmn, const PPM & ppm) {
        
        if(_gwa_energies.size()<1){
            throw std::runtime_error("Sigma gwa_energies not set!");
        }
      _sigma_x=Eigen::MatrixXd::Zero(_qptotal,_qptotal);
      _sigma_c=Eigen::MatrixXd::Zero(_qptotal,_qptotal);
      unsigned _levelsum = _Mmn.get_ntot(); // total number of bands
      unsigned _gwsize = _Mmn.getAuxDimension(); // size of the GW basis
      const double pi = boost::math::constants::pi<double>();

#pragma omp parallel for
      for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
        const MatrixXfd & Mmn = _Mmn[ _gw_level + _qpmin ];
        double sigma_x = 0;
        for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
          // loop over all occupied bands used in screening
          for (unsigned _i_occ = 0; _i_occ <= _homo; _i_occ++) {
            sigma_x -= Mmn(_i_gw, _i_occ) * Mmn(_i_gw, _i_occ);
          } // occupied bands
        } // gwbasis functions             
        _sigma_x(_gw_level, _gw_level) = (1.0 - _ScaHFX) * sigma_x;
      }
      // initial _qp_energies are dft energies
      Eigen::VectorXd _qp_old = _gwa_energies;
      // only diagonal elements except for in final iteration
      for (unsigned _g_iter = 0; _g_iter < _g_sc_max_iterations; _g_iter++) {
        // loop over all GW levels
//#pragma omp parallel for
        for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
          const MatrixXfd & Mmn = _Mmn[ _gw_level + _qpmin ];
          const double qpmin = _qp_old(_gw_level + _qpmin);
          double sigma_c = 0.0;
          // loop over all functions in GW basis
          for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
            // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
            if (ppm.getPpm_weight()(_i_gw) < 1.e-9) {
              continue;
            }
            const double ppm_freq = ppm.getPpm_freq()(_i_gw);
            const double fac = ppm.getPpm_weight()(_i_gw) * ppm_freq;
            // loop over all bands
            for (unsigned _i = 0; _i < _levelsum; _i++) {
              double occ = 1.0;
              if (_i > _homo) occ = -1.0; // sign for empty levels
              // energy denominator
              const double _denom = qpmin - _qp_old(_i) + occ*ppm_freq;
              double _stab = 1.0;
              if (std::abs(_denom) < 0.25) {
                _stab = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
              }
              const double _factor = 0.5 * fac * _stab / _denom; //Hartree
              // sigma_c diagonal elements
              sigma_c += _factor * Mmn(_i_gw, _i) * Mmn(_i_gw, _i);
            }// bands
          }// GW functions
          _sigma_c(_gw_level, _gw_level) = sigma_c;
          // update _qp_energies
          _gwa_energies(_gw_level + _qpmin) = (*_dftenergies)(_gw_level + _qpmin) + sigma_c + _sigma_x(_gw_level, _gw_level) - (*_vxc)(_gw_level, _gw_level);
        }// all bands
        Eigen::VectorXd diff = _qp_old - _gwa_energies;
        bool energies_converged = true;

        int state = 0;
        double diff_max = diff.cwiseAbs().maxCoeff(&state);
        if (diff_max > _g_sc_limit) {
          energies_converged = false;
        }

        if (tools::globals::verbose) {
          double _DFTgap = (*_dftenergies)(_homo + 1) - (*_dftenergies)(_homo);
          double _QPgap = _gwa_energies(_homo + 1) - _gwa_energies(_homo);
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " G_Iteration: " << _g_iter + 1 << " shift=" << _QPgap - _DFTgap << " E_diff max=" << diff_max << " StateNo:" << state << flush;
        }
        double alpha = 0.0;
        _gwa_energies = (1 - alpha) * _gwa_energies + alpha*_qp_old;

        if (energies_converged) {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converged after " << _g_iter + 1 << " G iterations." << flush;
          break;
        } else if (_g_iter == _g_sc_max_iterations - 1) {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " G-self-consistency cycle not converged after " << _g_sc_max_iterations << " iterations." << flush;
          break;

        } else {
          _qp_old = _gwa_energies;
        }
      } // iterations
      return;
    }

    void Sigma::X_offdiag(const TCMatrix_gwbse& _Mmn){
      unsigned _gwsize = _Mmn.getAuxDimension();
      #pragma omp parallel for schedule(dynamic)
      for (unsigned _gw_level1 = 0; _gw_level1 < _qptotal; _gw_level1++) {
        const MatrixXfd Mmn1 = _Mmn[ _gw_level1 + _qpmin ].block(0,0,_gwsize,_homo+1);
        for (unsigned _gw_level2 = _gw_level1+1; _gw_level2 < _qptotal; _gw_level2++) {
          const MatrixXfd & Mmn2 = _Mmn[ _gw_level2 + _qpmin ];
          double sigma_x = 0;
          for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
            // loop over all occupied bands used in screening
            for (unsigned _i_occ = 0; _i_occ <= _homo; _i_occ++) {
              sigma_x -= Mmn1(_i_gw, _i_occ) * Mmn2(_i_gw, _i_occ);
            } // occupied bands
          } // gwbasis functions
          _sigma_x(_gw_level1, _gw_level2) = (1.0 - _ScaHFX) * sigma_x;
          _sigma_x(_gw_level2, _gw_level1) = (1.0 - _ScaHFX) * sigma_x;
        }
      }
      return;
    }

    void Sigma::C_offdiag(const TCMatrix_gwbse& _Mmn, const PPM& ppm){
            
      #pragma omp parallel 
      {
        unsigned lumo=_homo+1;
        const unsigned _levelsum = _Mmn.get_ntot(); // total number of bands
        const unsigned _gwsize = _Mmn.getAuxDimension(); // size of the GW basis
        const double fourpi = 4*boost::math::constants::pi<double>();
        
        const Eigen::VectorXd ppm_weight=ppm.getPpm_weight();
        const Eigen::VectorXd ppm_freqs=ppm.getPpm_freq();
        #pragma omp for schedule(dynamic)
        for (unsigned _gw_level1 = 0; _gw_level1 < _qptotal; _gw_level1++) {
        const MatrixXfd Mmn1=_Mmn[ _gw_level1 + _qpmin ];
        for (unsigned _gw_level2 = _gw_level1+1; _gw_level2 < _qptotal; _gw_level2++) {
          const MatrixXfd Mmn1xMmn2=_Mmn[ _gw_level2 + _qpmin ].cwiseProduct(Mmn1);
          const Eigen::VectorXd gwa_energies=_gwa_energies;
          const double qpmin1 = gwa_energies(_gw_level1 + _qpmin);
          const double qpmin2 = gwa_energies(_gw_level2 + _qpmin);
          double sigma_c=0;
          for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
            // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
            if (ppm_weight(_i_gw) < 1.e-9) {
              continue;
            }
            const double ppm_freq= ppm_freqs(_i_gw);
            const double fac =0.25* ppm_weight(_i_gw) * ppm_freq;
            double sigma_loc=0.0;
            // loop over occ screening levels
            for (unsigned _i = 0; _i <lumo; _i++) {              
              const double gwa_energy = gwa_energies(_i)-ppm_freq;
              // energy denominator
              const double _denom1 = qpmin1 - gwa_energy;
              double _stab1 = 1.0;
              if (std::abs(_denom1) < 0.25) {
                  _stab1 = 0.5 * (1.0 - std::cos(fourpi * _denom1));
              }
              
              const double _denom2 = qpmin2 - gwa_energy;
              double _stab2 = 1.0;
              if (std::abs(_denom2) < 0.25) {
                  _stab2 = 0.5 * (1.0 - std::cos(fourpi * _denom2));
              }
              const double factor=_stab1/_denom1+_stab2 / _denom2; //Hartree}
              sigma_loc+=Mmn1xMmn2(_i_gw,_i)*factor;
            }
            // loop over unocc screening levels
            for (unsigned _i = lumo; _i < _levelsum; _i++) {              
              const double gwa_energy = gwa_energies(_i)+ppm_freq;
              // energy denominator
              const double _denom1 = qpmin1 - gwa_energy;
              double _stab1 = 1.0;
              if (std::abs(_denom1) < 0.25) {
                  _stab1 = 0.5 * (1.0 - std::cos(fourpi * _denom1));
              }
              
              const double _denom2 = qpmin2 - gwa_energy;
              double _stab2 = 1.0;
              if (std::abs(_denom2) < 0.25) {
                  _stab2 = 0.5 * (1.0 - std::cos(fourpi * _denom2));
              }
              const double factor=_stab1/_denom1+_stab2 / _denom2; //Hartree}
              sigma_loc+=Mmn1xMmn2(_i_gw,_i)*factor;
            }
            sigma_c+=sigma_loc*fac;
          }
          _sigma_c(_gw_level1, _gw_level2) = sigma_c;
          _sigma_c(_gw_level2, _gw_level1) = sigma_c;
        }// GW row             
      }//GW col
      }
      
      return;
    }

    void Sigma::CalcOffDiagElements(const TCMatrix_gwbse& _Mmn, const PPM & ppm) {
     
      X_offdiag(_Mmn);
      C_offdiag(_Mmn, ppm);
      return;
    }

  }
};
