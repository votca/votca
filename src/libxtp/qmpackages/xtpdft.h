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

#ifndef __VOTCA_XTP_XTPDFT_H
#define __VOTCA_XTP_XTPDFT_H


#include <votca/ctp/apolarsite.h>
#include <votca/xtp/qmpackage.h>
#include <votca/xtp/dftengine.h>

#include <string>



namespace votca {
    namespace xtp {

        /**
            \brief Wrapper for the internal XTP DFT engine


         */
        class XTPDFT : public QMPackage {
        public:

            std::string getPackageName() {
                return "xtp";
            }

            void Initialize(tools::Property &options);

            bool WriteInputFile( Orbitals& orbitals);

            bool Run(Orbitals& orbitals);

            void CleanUp();

            bool CheckLogFile();

            bool ParseLogFile(Orbitals& orbitals);

            bool ParseOrbitalsFile(Orbitals& orbitals);
            
            void setMultipoleBackground( std::vector<std::shared_ptr<ctp::PolarSeg> > multipoles);

        private:
            void WriteChargeOption() { return ;}
            tools::Property _xtpdft_options;
            std::string _cleanup;

            
        };


    }
}

#endif /* __VOTCA_XTP_XTPDFT_H */
