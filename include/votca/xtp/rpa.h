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

#ifndef _VOTCA_XTP_RPA_H
#define _VOTCA_XTP_RPA_H
#include <votca/xtp/votca_config.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/threecenter.h>

namespace votca {
    namespace xtp {

        class RPA {
        public:

            void configure(unsigned homo, unsigned rpamin, unsigned rpamax) {
                _homo = homo;
                _rpamin = rpamin;
                _rpamax = rpamax;
                // fix the frequencies for PPM
                _screening_freq = Eigen::MatrixXd::Zero(2, 2); // two frequencies
                // first one
                _screening_freq(0, 0) = 0.0; // real part
                _screening_freq(0, 1) = 0.0; // imaginary part
                // second one
                _screening_freq(1, 0) = 0.0; // real part
                _screening_freq(1, 1) = 0.5; // imaginary part  //hartree
                // one entry to epsilon for each frequency
                _epsilon.resize(_screening_freq.rows());
            }

            const Eigen::MatrixXd& GetScreening_freq() const {
                return _screening_freq;
            }

            const std::vector<Eigen::MatrixXd>& GetEpsilon() const {
                return _epsilon;
            }

            void prepare_threecenters(const TCMatrix_gwbse& _Mmn_full);

            void calculate_epsilon(const Eigen::VectorXd& qp_energies);

            void FreeMatrices() {
                _Mmn_RPA.Cleanup();
                _epsilon[0].resize(0, 0);
                _epsilon[1].resize(0, 0);

            }

        private:
            TCMatrix_gwbse _Mmn_RPA;

            unsigned _homo; // HOMO index
            unsigned _rpamin;
            unsigned _rpamax;
            double _shift; // pre-shift of DFT energies

            // container for the epsilon matrix
            std::vector<Eigen::MatrixXd > _epsilon;
            // container for frequencies in screening (index 0: real part, index 1:
            // imaginary part)
            Eigen::MatrixXd _screening_freq;

            Eigen::MatrixXd epsilon_real(const Eigen::VectorXd& qp_energies,
                    const double screening_freq);

            Eigen::MatrixXd epsilon_imaginary(const Eigen::VectorXd& qp_energies,
                    const double screening_freq);


        };
    }
}

#endif /* _VOTCA_RPA_RPA_H */
