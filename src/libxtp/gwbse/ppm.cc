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

#include "votca/xtp/ppm.h"



namespace votca {
    namespace xtp {


        void PPM::PPM_construct_parameters(const RPA& rpa,const Eigen::MatrixXd& _overlap_cholesky_inverse) {
            
            // orthogonalize via L-1 epsilon L-T
             
            Eigen::MatrixXd ortho = _overlap_cholesky_inverse*rpa.GetEpsilon()[0]*_overlap_cholesky_inverse.transpose();
            //Solve Eigensystem
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho); 
            //we store _ppm_phi_T instead of _ppm_phi because we need it for later transformations
            _ppm_phi_T=es.eigenvectors().transpose()*_overlap_cholesky_inverse;
            
            // store PPM weights from eigenvalues
            _ppm_weight.resize(es.eigenvalues().size());
            for (unsigned _i = 0; _i < es.eigenvalues().size(); _i++) {
                _ppm_weight(_i) = 1.0 - 1.0 / es.eigenvalues()(_i);
            }

            // determine PPM frequencies
            _ppm_freq.resize(es.eigenvalues().size());
            // a) phi^t * epsilon(1) * phi e.g. transform epsilon(1) to the same space as epsilon(0)
           ortho=_ppm_phi_T*_epsilon[1]*_ppm_phi_T.transpose();
           Eigen::MatrixXd epsilon_1_inv=ortho.inverse();
           
            
            #pragma omp parallel for 
            for (unsigned _i = 0; _i < es.eigenvalues().size(); _i++) {

                if (rpa.GetScreening_freq()(1, 0) == 0.0) {
                    if (_ppm_weight(_i) < 1.e-5) {
                        _ppm_weight(_i) = 0.0;
                        _ppm_freq(_i) = 0.5;//Hartree
                        continue;
                    } else {
                        double _nom = epsilon_1_inv(_i, _i) - 1.0;
                        double _frac = -1.0 * _nom / (_nom + _ppm_weight(_i)) * rpa.GetScreening_freq()(1, 1) * rpa.GetScreening_freq()(1, 1);
                        _ppm_freq(_i) = sqrt(std::abs(_frac));
                    }

                } else {
                    // only purely imaginary frequency assumed
                    cerr << " mixed frequency! real part: " << rpa.GetScreening_freq()(1, 0) << " imaginary part: " << rpa.GetScreening_freq()(1, 1) << flush;
                    exit(1);
                }

            }
 
           
            

            return;
        }

 
}};
