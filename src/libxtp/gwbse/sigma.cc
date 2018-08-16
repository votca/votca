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

    void Sigma::X_diag(const TCMatrix_gwbse& Mmn){
      unsigned levelsum = Mmn.get_ntot(); // total number of bands
      unsigned gwsize = Mmn.getAuxDimension(); // size of the GW basis
      _sigma_x=Eigen::MatrixXd::Zero(_qptotal,_qptotal);
      #pragma omp parallel for
      for (unsigned gw_level = 0; gw_level < _qptotal; gw_level++) {
        const MatrixXfd & Mmn1 = Mmn[ gw_level + _qpmin ];
        double sigma_x = 0;
        for (unsigned i_gw = 0; i_gw < gwsize; i_gw++) {
          // loop over all occupied bands used in screening
          for (unsigned i_occ = 0; i_occ <= _homo; i_occ++) {
            sigma_x -= Mmn1( i_occ,i_gw) * Mmn1( i_occ,i_gw);
          } // occupied bands
        } // gwbasis functions             
        _sigma_x(gw_level, gw_level) = (1.0 - _ScaHFX) * sigma_x;
      }
    }

    void Sigma::C_diag(const TCMatrix_gwbse& Mmn, const PPM& ppm, const Eigen::VectorXd& qp_old){
      unsigned levelsum = Mmn.get_ntot(); // total number of bands
      unsigned gwsize = Mmn.getAuxDimension(); // size of the GW basis
      const double fourpi = 4*boost::math::constants::pi<double>();
      
      // loop over all GW levels
#pragma omp parallel for
      for (unsigned gw_level = 0; gw_level < _qptotal; gw_level++) {
        const MatrixXfd & Mmn1 = Mmn[ gw_level + _qpmin ];
        const double qpmin = qp_old(gw_level + _qpmin);
        double sigma_c = 0.0;
        // loop over all functions in GW basis
        for (unsigned i_gw = 0; i_gw < gwsize; i_gw++) {
          // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
          if (ppm.getPpm_weight()(i_gw) < 1.e-9) {
            continue;
          }
          const double ppm_freq = ppm.getPpm_freq()(i_gw);
          const double fac = 0.5*ppm.getPpm_weight()(i_gw) * ppm_freq;
          // loop over all bands
          double sigma_c_loc=0.0;
          for (unsigned i = 0; i < _homo+1; i++) {
            const double denom = qpmin - qp_old(i) +ppm_freq;
            double stab = 1.0;
            if (std::abs(denom) < 0.25) {
              stab = 0.5 * (1.0 - std::cos(fourpi * denom));
            }
            const double _factor = stab / denom; //Hartree
            // sigma_c diagonal elements
            sigma_c_loc += _factor * Mmn1(i, i_gw) * Mmn1( i,i_gw);
          }// bands
          for (unsigned i = _homo+1; i < levelsum; i++) {
            const double denom = qpmin - qp_old(i) -ppm_freq;
            double stab = 1.0;
            if (std::abs(denom) < 0.25) {
              stab = 0.5 * (1.0 - std::cos(fourpi * denom));
            }
            const double factor = stab / denom; //Hartree
            // sigma_c diagonal elements
            sigma_c_loc += factor * Mmn1(i,i_gw) * Mmn1(i,i_gw);
          }// bands
          sigma_c+=sigma_c_loc*fac;
        }// GW functions
        _sigma_c(gw_level, gw_level) = sigma_c;
        // update _qp_energies
        _gwa_energies(gw_level + _qpmin) = (*_dftenergies)(gw_level + _qpmin) + sigma_c + _sigma_x(gw_level, gw_level) - (*_vxc)(gw_level, gw_level);
      }// all bands
    }

    void Sigma::CalcdiagElements(const TCMatrix_gwbse& Mmn, const PPM & ppm) {
        X_diag(Mmn);
        if(_gwa_energies.size()<1){
            throw std::runtime_error("Sigma gwa_energies not set!");
        }
      _sigma_c=Eigen::MatrixXd::Zero(_qptotal,_qptotal);
      


      // initial _qp_energies are dft energies
      Eigen::VectorXd qp_old = _gwa_energies;
      // only diagonal elements except for in final iteration
      for (unsigned g_iter = 0; g_iter < _g_sc_max_iterations; g_iter++) {
        
        C_diag(Mmn, ppm, qp_old);
        Eigen::VectorXd diff = qp_old - _gwa_energies;
        bool energies_converged = true;

        int state = 0;
        double diff_max = diff.cwiseAbs().maxCoeff(&state);
        if (diff_max > _g_sc_limit) {
          energies_converged = false;
        }

        if (tools::globals::verbose) {
          double _DFTgap = (*_dftenergies)(_homo + 1) - (*_dftenergies)(_homo);
          double _QPgap = _gwa_energies(_homo + 1) - _gwa_energies(_homo);
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " G_Iteration: " << g_iter + 1 << " shift=" << _QPgap - _DFTgap << " E_diff max=" << diff_max << " StateNo:" << state << std::flush;
        }
        double alpha = 0.0;
        _gwa_energies = (1 - alpha) * _gwa_energies + alpha*qp_old;

        if (energies_converged) {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converged after " << g_iter + 1 << " G iterations." << std::flush;
          break;
        } else if (g_iter == _g_sc_max_iterations - 1) {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " G-self-consistency cycle not converged after " << _g_sc_max_iterations << " iterations." << std::flush;
          break;

        } else {
          qp_old = _gwa_energies;
        }
      } // iterations
      return;
    }

    void Sigma::X_offdiag(const TCMatrix_gwbse& Mmn){
      unsigned _gwsize = Mmn.getAuxDimension();
      #pragma omp parallel for schedule(dynamic)
      for (unsigned gw_level1 = 0; gw_level1 < _qptotal; gw_level1++) {
        const MatrixXfd& Mmn1 = Mmn[ gw_level1 + _qpmin ];
        for (unsigned gw_level2 = gw_level1+1; gw_level2 < _qptotal; gw_level2++) {
          const MatrixXfd & Mmn2 = Mmn[ gw_level2 + _qpmin ];
          double sigma_x = 0;
          for (unsigned i_gw = 0; i_gw < _gwsize; i_gw++) {
            // loop over all occupied bands used in screening
            for (unsigned i_occ = 0; i_occ <= _homo; i_occ++) {
              sigma_x -= Mmn1(i_occ,i_gw) * Mmn2(i_occ,i_gw);
            } // occupied bands
          } // gwbasis functions
          _sigma_x(gw_level1, gw_level2) = (1.0 - _ScaHFX) * sigma_x;
          _sigma_x(gw_level2, gw_level1) = (1.0 - _ScaHFX) * sigma_x;
        }
      }
      return;
    }

    double Sigma::SumSymmetric(real_gwbse Mmn1xMmn2, double qpmin1, double qpmin2, double gwa_energy){
      const double fourpi = 4*boost::math::constants::pi<double>();
      double denom = qpmin1 - gwa_energy;
      double stab = 1.0;
      if (std::abs(denom) < 0.25) {
        stab = 0.5 * (1.0 - std::cos(fourpi * denom));
      }
      double factor= stab / denom;
      denom = qpmin2 - gwa_energy;
      stab = 1.0;
      if (std::abs(denom) < 0.25) {
        stab = 0.5 * (1.0 - std::cos(fourpi * denom));
      }
      factor+= stab / denom; //Hartree}
      return Mmn1xMmn2 * factor;
    }

    void Sigma::C_offdiag(const TCMatrix_gwbse& Mmn, const PPM& ppm){
            
      #pragma omp parallel 
      {
        unsigned lumo=_homo+1;
        const unsigned levelsum = Mmn.get_ntot(); // total number of bands
        const unsigned gwsize = Mmn.getAuxDimension(); // size of the GW basis
        const Eigen::VectorXd ppm_weight=ppm.getPpm_weight();
        const Eigen::VectorXd ppm_freqs=ppm.getPpm_freq();
        #pragma omp for schedule(dynamic)
        for (unsigned gw_level1 = 0; gw_level1 < _qptotal; gw_level1++) {
        const MatrixXfd& Mmn1=Mmn[ gw_level1 + _qpmin ];
        for (unsigned gw_level2 = gw_level1+1; gw_level2 < _qptotal; gw_level2++) {
          const MatrixXfd Mmn1xMmn2=Mmn[ gw_level2 + _qpmin ].cwiseProduct(Mmn1);
          const Eigen::VectorXd gwa_energies=_gwa_energies;
          const double qpmin1 = gwa_energies(gw_level1 + _qpmin);
          const double qpmin2 = gwa_energies(gw_level2 + _qpmin);
 
          double sigma_c=0;
          for (unsigned i_gw = 0; i_gw < gwsize; i_gw++) {
            // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
            if (ppm_weight(i_gw) < 1.e-9) {
              continue;
            }
            const double ppm_freq= ppm_freqs(i_gw);
            const double fac = 0.25* ppm_weight(i_gw) * ppm_freq;
            double sigma_loc=0.0;
            // loop over occ screening levels
              for (unsigned i = 0; i < lumo; i++) {
                const double gwa_energy = gwa_energies(i) - ppm_freq;
                sigma_loc+=SumSymmetric(Mmn1xMmn2(i,i_gw), qpmin1, qpmin2, gwa_energy);
              }
              // loop over unocc screening levels
              for (unsigned i = lumo; i < levelsum; i++) {
                const double gwa_energy = gwa_energies(i) + ppm_freq;
                sigma_loc+=SumSymmetric(Mmn1xMmn2(i,i_gw), qpmin1, qpmin2, gwa_energy);
              }
              sigma_c += sigma_loc*fac;
            }
          _sigma_c(gw_level1, gw_level2) = sigma_c;
          _sigma_c(gw_level2, gw_level1) = sigma_c;
        }// GW row             
      }//GW col
      }
      
      return;
    }

    void Sigma::CalcOffDiagElements(const TCMatrix_gwbse& Mmn, const PPM & ppm) {
     
      X_offdiag(Mmn);
      C_offdiag(Mmn, ppm);

      return;
    }

  }
};
