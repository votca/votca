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

#ifndef VOTCA_XTP_DFTENGINE_H
#define VOTCA_XTP_DFTENGINE_H

#include <votca/xtp/numerical_integrations.h>
#include <boost/filesystem.hpp>
#include <votca/xtp/logger.h>
#include <votca/xtp/polarsegment.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/ERIs.h>
#include <votca/xtp/convergenceacc.h>

namespace votca {
    namespace xtp {
        class Orbitals;

/**
         * \brief Electronic ground-state via Density-Functional Theory
         *
         * Evaluates electronic ground state in molecular systems based on
         * density functional theory with Gaussian Orbitals.
         * 
         */

        class DFTEngine {
        public:

            DFTEngine(Orbitals& orbitals):_orbitals(orbitals){};

         
            void Initialize(tools::Property &options);

            void CleanUp();

            void setLogger(Logger* pLog) {
                _pLog = pLog;
            }

            void ConfigureExternalGrid(const std::string& grid_name_ext) {
                _grid_name_ext = grid_name_ext;
                _do_externalfield = true;
            }

            void setExternalcharges(std::shared_ptr<MMRegion>& externalsites) {
                _externalsites = externalsites;           
                _addexternalsites = true;
            }
                   
            void setExternalGrid(std::vector<double> electrongrid, std::vector<double> nucleigrid) {
                _externalgrid = electrongrid;
                _externalgrid_nuc = nucleigrid;
            }

            std::vector< const Eigen::Vector3d *> getExternalGridpoints() {
                return _gridIntegration_ext.getGridpoints();
            }

            bool Evaluate();
   
            void Prepare();

            std::string getDFTBasisName() const{
                return _dftbasis_name;
            };
            
            

        private:
            
            Eigen::MatrixXd OrthogonalizeGuess(const Eigen::MatrixXd& GuessMOs )const;
            void PrintMOs(const Eigen::VectorXd& MOEnergies);
            void CalcElDipole()const;
            void CalculateERIs(const AOBasis& dftbasis, const Eigen::MatrixXd &DMAT);
            void ConfigOrbfile();
            void SetupInvariantMatrices();
            Eigen::MatrixXd AtomicGuess();
            std::string ReturnSmallGrid(const std::string& largegrid);
            
            Eigen::MatrixXd IntegrateExternalDensity(const Orbitals& extdensity);
            
            Eigen::MatrixXd RunAtomicDFT_unrestricted(const QMAtom& uniqueAtom);
            
            void NuclearRepulsion();
            double ExternalRepulsion();
            double ExternalGridRepulsion(std::vector<double> externalpotential_nuc);
            Eigen::MatrixXd SphericalAverageShells(const Eigen::MatrixXd& dmat, AOBasis& dftbasis);

            Logger *_pLog;

            int _openmp_threads;

            // atoms
            Orbitals& _orbitals;

            // basis sets
            std::string _auxbasis_name;
            std::string _dftbasis_name;
            std::string _ecp_name;
            BasisSet _dftbasisset;
            BasisSet _auxbasisset;
            BasisSet _ecpbasisset;
            AOBasis _dftbasis;
            AOBasis _auxbasis;
            AOBasis _ecp;

            bool _with_ecp;
            bool _with_RI;
            
            std::string _four_center_method; // direct | cache
            
            // Pre-screening
            bool _with_screening;
            double _screening_eps;
            
            // numerical integration Vxc
            std::string _grid_name;
            std::string _grid_name_small;
            bool _use_small_grid;
            NumericalIntegration _gridIntegration;
            NumericalIntegration _gridIntegration_small;
            //used to store Vxc after final iteration

            //numerical integration externalfield;
            //this will not remain here but be moved to qmape
            bool _do_externalfield=false;
            std::string _grid_name_ext;
            NumericalIntegration _gridIntegration_ext;
            std::vector<double> _externalgrid;
            std::vector<double> _externalgrid_nuc;

            Eigen::MatrixXd _dftAOdmat;

            // AO Matrices
            AOOverlap _dftAOoverlap;
            AOKinetic _dftAOkinetic;
            AOESP _dftAOESP;
            AOECP _dftAOECP;
            AODipole_Potential _dftAODipole_Potential;
            AOQuadrupole_Potential _dftAOQuadrupole_Potential;
            AOPlanewave _dftAOplanewave;
            double _E_nucnuc;
            
            bool _with_guess;
            std::string _initial_guess;

            // Convergence
            int _numofelectrons=0;
            int _max_iter=100;
            ConvergenceAcc::options _conv_opt;
            //DIIS variables
            ConvergenceAcc _conv_accelerator;
            //Electron repulsion integrals
            ERIs _ERIs;

            // external charges
            std::shared_ptr<MMRegion> _externalsites;
            bool _addexternalsites=false;

            // exchange and correlation
            double _ScaHFX;
            std::string _xc_functional_name;

            bool _integrate_ext_density=false;
            //integrate external density
            std::string _orbfilename;
            std::string _gridquality;
            std::string _state;
        };


    }
}

#endif // VOTCA_XTP_DFTENGINE_H
