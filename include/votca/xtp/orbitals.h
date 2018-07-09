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

#ifndef __VOTCA_XTP_ORBITALS_H
#define __VOTCA_XTP_ORBITALS_H


//for openmp
#include <votca/xtp/eigen.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/qmatom.h>

#include <votca/xtp/checkpoint.h>
#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

#include <votca/ctp/logger.h>
#include <boost/format.hpp>
#include <votca/tools/constants.h>



namespace votca {
    namespace xtp {
       

        /**
         * \brief container for molecular orbitals
         *
         * The Orbitals class stores orbital id, energy, MO coefficients, basis set
         *
         */
        class Orbitals {
        public:

            Orbitals();
            ~Orbitals();

            static void PrepareGuess(Orbitals* _orbitalsA, Orbitals* _orbitalsB, Orbitals* _orbitalsAB);

            /*
             *
             *  ************** NEW ACCESS STRATEGY ****************
             *
             *  Scalars:              get and set functions
             *  Vectors and matrixes: const and non-const refs, has-function via size
             */

            // access to DFT basis set size, new, tested

            bool hasBasisSetSize() {
                return ( _basis_set_size > 0) ? true : false;
            }

            int getBasisSetSize() {
                return _basis_set_size;
            }

            void setBasisSetSize(const int &basis_set_size) {
                _basis_set_size = basis_set_size;
            }

            // access to DFT number of levels, new, tested

            bool hasNumberOfLevels() {
                return ( (_occupied_levels > 0) && (_unoccupied_levels > 0) ? true : false);
            }

            unsigned int getNumberOfLevels() {
                return ( _occupied_levels + _unoccupied_levels);
            }
            void setNumberOfLevels(const int &occupied_levels, const int &unoccupied_levels);

            // access to DFT number of electrons, new, tested

            bool hasNumberOfElectrons() {
                return ( _number_of_electrons > 0) ? true : false;
            }

            int getNumberOfElectrons() const{
                return _number_of_electrons;
            };

            void setNumberOfElectrons(int electrons) {
                _number_of_electrons = electrons;
            }


            bool hasECP(){
                return ( _ECP !="") ? true : false;
            }

            std::string getECP() {
                return _ECP;
            };

            void setECP(const std::string &ECP) {
                _ECP = ECP;
            };

            // access to QM package name, new, tested

            bool hasQMpackage() {
                return (!_qm_package.empty());
            }

            std::string getQMpackage() {
                return _qm_package;
            }

            void setQMpackage(std::string qmpackage) {
                _qm_package = qmpackage;
            }

            // access to DFT AO overlap matrix, new, tested

            bool hasAOOverlap() {
                return ( _overlap.rows() > 0) ? true : false;
            }

            const Eigen::MatrixXd &AOOverlap() const {
                return _overlap;
            }

            Eigen::MatrixXd &AOOverlap() {
                return _overlap;
            }

            // access to DFT molecular orbital energies, new, tested

            bool hasMOEnergies() {
                return ( _mo_energies.size() > 0) ? true : false;
            }

            const Eigen::VectorXd &MOEnergies() const {
                return _mo_energies;
            }

            Eigen::VectorXd &MOEnergies() {
                return _mo_energies;
            }

            // access to DFT molecular orbital energy of a specific level (in eV)

            double getEnergy(int level) {
                return ( hasMOEnergies()) ? votca::tools::conv::hrt2ev * _mo_energies[level - 1] : 0;
            }

            // access to DFT molecular orbital coefficients, new, tested

            bool hasMOCoefficients() {
                return ( _mo_coefficients.cols() > 0) ? true : false;
            }

            const Eigen::MatrixXd &MOCoefficients() const {
                return _mo_coefficients;
            }

            Eigen::MatrixXd &MOCoefficients() {
                return _mo_coefficients;
            }

            // access to DFT transfer integrals, new, tested

            bool hasMOCouplings() {
                return ( _mo_couplings.cols() > 0) ? true : false;
            }

            const Eigen::MatrixXd &MOCouplings() const {
                return _mo_couplings;
            }

            Eigen::MatrixXd &MOCouplings() {
                return _mo_couplings;
            }


            // determine (pseudo-)degeneracy of a DFT molecular orbital

            bool hasDegeneracy() {
                return ( !_level_degeneracy.empty()) ? true : false;
            }
            std::vector<int>* getDegeneracy(int level, double _energy_difference);

            // access to QM atoms

            bool hasQMAtoms() {
                return ( _atoms.size() > 0) ? true : false;
            }

            const std::vector< QMAtom* > &QMAtoms() const {
                return _atoms;
            }

            std::vector< QMAtom* > &QMAtoms() {
                return _atoms;
            }

            // access to classical self-energy in MM environment, new, tested

            bool hasSelfEnergy() {
                return ( _self_energy != 0.0) ? true : false;
            }

            double getSelfEnergy() {
                return _self_energy;
            }

            void setSelfEnergy(double selfenergy) {
                _self_energy = selfenergy;
            }

            // access to QM total energy, new, tested

            bool hasQMEnergy() {
                return ( _qm_energy != 0.0) ? true : false;
            }

            double getQMEnergy() {
                return _qm_energy;
            }

            void setQMEnergy(double qmenergy) {
                _qm_energy = qmenergy;
            }

            // access to DFT basis set name

            bool hasDFTbasis() {
                return ( !_dftbasis.empty()) ? true : false;
            }

            void setDFTbasis(const std::string basis) {
                _dftbasis = basis;
            }

            const std::string getDFTbasis() const {
                return _dftbasis;
            }



            /*
             *  ======= GW-BSE related functions =======
             */

            // access to exchange-correlation AO matrix, new, tested

            bool hasAOVxc() {
                return ( _vxc.rows() > 0) ? true : false;
            }

            Eigen::MatrixXd &AOVxc() {
                return _vxc;
            }

            const Eigen::MatrixXd &AOVxc() const {
                return _vxc;
            }

            // access to auxiliary basis set name

            bool hasAuxbasis() {
                return ( !_auxbasis.empty()) ? true : false;
            }

            void setAuxbasis(std::string basis) {
                _auxbasis = basis;
            }

            const std::string getAuxbasis() const {
                return _auxbasis;
            }


            // access to list of indices used in GWA

            bool hasGWAindices() {
                return ( _qpmax > 0) ? true : false;
            }

            void setGWAindices(int qpmin, int qpmax) {
                _qpmin = qpmin;
                _qpmax = qpmax;
                _qptotal = _qpmax - _qpmin + 1;
            }

            unsigned getGWAmin() const {
                return _qpmin;
            }

            unsigned getGWAmax() const {
                return _qpmax;
            }

            unsigned getGWAtot() const {
                return (_qpmax - _qpmin + 1);
            }

            // access to list of indices used in RPA

            bool hasRPAindices() {
                return ( _rpamax > 0) ? true : false;
            }

            void setRPAindices(int rpamin, int rpamax) {
                _rpamin = rpamin;
                _rpamax = rpamax;
            }

            int getRPAmin() const {
                return _rpamin;
            }

            int getRPAmax() const {
                return _rpamax;
            }

            // access to list of indices used in BSE

            void setBSEtype(std::string bsetype){_bsetype=bsetype;}
            std::string getBSEtype() const{return _bsetype;}


            bool hasBSEindices() {
                return ( _bse_cmax > 0) ? true : false;
            }

            void setBSEindices(int vmin, int vmax, int cmin, int cmax, int nmax) {
                _bse_vmin = vmin;
                _bse_vmax = vmax;
                _bse_cmin = cmin;
                _bse_cmax = cmax;
                _bse_nmax = nmax;
                _bse_vtotal = _bse_vmax - _bse_vmin + 1;
                _bse_ctotal = _bse_cmax - _bse_cmin + 1;
                _bse_size = _bse_vtotal * _bse_ctotal;
                return;
            }

            int getBSEvmin() const {
                return _bse_vmin;
            }

            int getBSEvmax() const {
                return _bse_vmax;
            }

            int getBSEcmin() const {
                return _bse_cmin;
            }

            int getBSEcmax() const {
                return _bse_cmax;
            }

            double getScaHFX() const {
                return _ScaHFX;
            }

            void setScaHFX(double ScaHFX) {
                _ScaHFX = ScaHFX;
            }

            // access to perturbative QP energies

            bool hasQPpert() {
                return ( _QPpert_energies.size() > 0) ? true : false;
            }

            const Eigen::MatrixXd &QPpertEnergies() const {
                return _QPpert_energies;
            }

            Eigen::MatrixXd &QPpertEnergies() {
                return _QPpert_energies;
            }

            // access to diagonalized QP energies and wavefunctions

            bool hasQPdiag() {
                return ( _QPdiag_energies.size() > 0) ? true : false;
            }

            const Eigen::VectorXd &QPdiagEnergies() const {
                return _QPdiag_energies;
            }

            Eigen::VectorXd &QPdiagEnergies() {
                return _QPdiag_energies;
            }

            const Eigen::MatrixXd &QPdiagCoefficients() const {
                return _QPdiag_coefficients;
            }

            Eigen::MatrixXd &QPdiagCoefficients() {
                return _QPdiag_coefficients;
            }

            // access to eh interaction



            bool hasEHinteraction_triplet() {
                return ( _eh_t.cols() > 0) ? true : false;
            }
            
            bool hasEHinteraction_singlet() {
                return ( _eh_s.cols() > 0) ? true : false;
            }

            const MatrixXfd &eh_s() const {
                return _eh_s;
            }

            MatrixXfd &eh_s() {
                return _eh_s;
            }

            const MatrixXfd &eh_t() const {
                return _eh_t;
            }

            MatrixXfd &eh_t() {
                return _eh_t;
            }

            // access to triplet energies and wave function coefficients

            bool hasBSETriplets() {
                return ( _BSE_triplet_energies.cols() > 0) ? true : false;
            }

            const VectorXfd &BSETripletEnergies() const {
                return _BSE_triplet_energies;
            }

            VectorXfd &BSETripletEnergies() {
                return _BSE_triplet_energies;
            }

            const MatrixXfd &BSETripletCoefficients() const {
                return _BSE_triplet_coefficients;
            }

            MatrixXfd &BSETripletCoefficients() {
                return _BSE_triplet_coefficients;
            }

            // access to singlet energies and wave function coefficients

            bool hasBSESinglets() {
                return (_BSE_singlet_energies.cols() > 0) ? true : false;
            }

            const VectorXfd &BSESingletEnergies() const {
                return _BSE_singlet_energies;
            }

            VectorXfd &BSESingletEnergies() {
                return _BSE_singlet_energies;
            }

            const MatrixXfd &BSESingletCoefficients() const {
                return _BSE_singlet_coefficients;
            }

            MatrixXfd &BSESingletCoefficients() {
                return _BSE_singlet_coefficients;
            }

            // for anti-resonant part in full BSE

            const MatrixXfd &BSESingletCoefficientsAR() const {
                return _BSE_singlet_coefficients_AR;
            }

            MatrixXfd &BSESingletCoefficientsAR() {
                return _BSE_singlet_coefficients_AR;
            }

            // access to transition dipole moments

            bool hasTransitionDipoles() {
                return (_transition_dipoles.size() > 0) ? true : false;
            }

            const std::vector< tools::vec > &TransitionDipoles() const {
                return _transition_dipoles;
            }

            std::vector< tools::vec > &TransitionDipoles() {
                return _transition_dipoles;
            }

            std::vector<double> Oscillatorstrengths();

            // access to singlet coupling elements

            bool hasSingletCouplings() {
                return (_BSE_singlet_couplings.cols() > 0) ? true : false;
            }

            const Eigen::MatrixXd &SingletCouplings() const {
                return _BSE_singlet_couplings;
            }

            Eigen::MatrixXd &SingletCouplings() {
                return _BSE_singlet_couplings;
            }

            void setSingletCouplings(Eigen::MatrixXd couplings) {
                _BSE_singlet_couplings = couplings;
            }

            // access to triplet coupling elements

            bool hasTripletCouplings() {
                return (_BSE_triplet_couplings.cols() > 0) ? true : false;
            }

            const Eigen::MatrixXd &TripletCouplings() const {
                return _BSE_triplet_couplings;
            }

            Eigen::MatrixXd &TripletCouplings() {
                return _BSE_triplet_couplings;
            }

            void setTripletCouplings(Eigen::MatrixXd couplings) {
                _BSE_triplet_couplings = couplings;
            }

            // exciton coupling number of levels information

            bool hasCoupledExcitonsA() {
                return ( _couplingsA > 0) ? true : false;
            }

            int getCoupledExcitonsA() {
                return _couplingsA;
            }

            void setCoupledExcitonsA(const int &excitons) {
                _couplingsA = excitons;
            }

            bool hasCoupledExcitonsB() {
                return ( _couplingsB > 0) ? true : false;
            }

            int getCoupledExcitonsB() {
                return _couplingsB;
            }

            void setCoupledExcitonsB(const int &excitons) {
                _couplingsB = excitons;
            }


            // functions for calculating density matrices
            Eigen::MatrixXd DensityMatrixGroundState();
            std::vector<Eigen::MatrixXd > DensityMatrixExcitedState(const std::string& spin,int state = 0);
            Eigen::MatrixXd TransitionDensityMatrix(const std::string& spin,int state = 0);
            Eigen::MatrixXd DensityMatrixQuasiParticle(int state = 0);
            Eigen::MatrixXd LambdaMatrixQuasiParticle();



            double GetTotalEnergy(std::string _spintype, int _opt_state);

            // functions for analyzing fragment charges via Mulliken populations
            Eigen::VectorXd LoewdinPopulation(const Eigen::MatrixXd& _densitymatrix, const Eigen::MatrixXd& _overlapmatrix, int _frag);

            // access to fragment charges of singlet excitations

            bool hasFragmentChargesSingEXC() {
                return (_DqS_frag.size() > 0) ? true : false;
            }

            const std::vector< Eigen::VectorXd > &getFragmentChargesSingEXC() const {
                return _DqS_frag;
            }

             void setFragmentChargesSingEXC(std::vector< Eigen::VectorXd > DqS_frag) {
                _DqS_frag=DqS_frag;
            }



            // access to fragment charges of triplet excitations

            bool hasFragmentChargesTripEXC() {
                return (_DqT_frag.size() > 0) ? true : false;
            }

            const std::vector< Eigen::VectorXd > &getFragmentChargesTripEXC() const {
                return _DqT_frag;
            }

            void setFragmentChargesTripEXC(std::vector< Eigen::VectorXd > DqT_frag) {
                _DqT_frag=DqT_frag;
            }

            // access to fragment charges in ground state

            const Eigen::VectorXd &getFragmentChargesGS() const {
                return _GSq_frag;
            }

             void setFragmentChargesGS(Eigen::VectorXd GSq_frag) {
                 _GSq_frag=GSq_frag;
            }

            void setFragment_E_localisation_singlet(std::vector< Eigen::VectorXd >& popE){
                _popE_s=popE;
            }

            void setFragment_H_localisation_singlet(std::vector< Eigen::VectorXd > & popH){
                _popH_s=popH;
            }

            void setFragment_E_localisation_triplet(std::vector< Eigen::VectorXd > & popE){
                _popE_t=popE;
            }

            void setFragment_H_localisation_triplet(std::vector< Eigen::VectorXd > & popH){
                _popE_s=popH;
            }


            const std::vector< Eigen::VectorXd >& getFragment_E_localisation_singlet()const{
                return _popE_s;
            }
            const std::vector< Eigen::VectorXd >& getFragment_H_localisation_singlet()const{
                return _popH_s;
            }
            const std::vector< Eigen::VectorXd >& getFragment_E_localisation_triplet()const{
                return _popE_t;
            }
            const std::vector< Eigen::VectorXd >& getFragment_H_localisation_triplet()const{
                return _popH_t;
            }

            Eigen::VectorXd FragmentNuclearCharges(int _frag);

            // returns indeces of a re-sorted in a descending order vector of energies
            std::vector<int> SortEnergies();

            /** Adds a QM atom to the atom list */
            QMAtom* AddAtom(int AtomID,std::string _type,double _x, double _y, double _z,
                    double _charge = 0) {
                QMAtom* pAtom = new QMAtom(AtomID,_type, _x, _y, _z);
                _atoms.push_back(pAtom);
                return pAtom;
            }

            QMAtom* AddAtom(int AtomID,std::string _type, tools::vec pos) {
                QMAtom* pAtom = new QMAtom(AtomID,_type, pos);
                _atoms.push_back(pAtom);
                return pAtom;
            }

            QMAtom* AddAtom(QMAtom atom) {
                QMAtom* pAtom = new QMAtom(atom);
                _atoms.push_back(pAtom);
                return pAtom;
            }


            void WritePDB(FILE *out, std::string tag = "");

            // reduces number of virtual orbitals to factor*number_of_occupied_orbitals
            void Trim(int factor);

            // reduces number of virtual orbitals to [HOMO-degG:LUMO+degL]
            void Trim(int degH, int degL);

            void LoadFromXYZ(std::string filename);

            void WriteToCpt(const std::string& filename);
            
            void ReadFromCpt(const std::string& filename);
            


        private:
            
            void WriteToCpt(CheckpointFile f);
            void WriteToCpt(CptLoc parent);
            
            void ReadFromCpt(CheckpointFile f);
            void ReadFromCpt(CptLoc parent);
            std::vector<Eigen::MatrixXd > DensityMatrixExcitedState_R(const std::string& spin,int state = 0);
            std::vector<Eigen::MatrixXd >DensityMatrixExcitedState_AR(const std::string& spin,int state = 0);

            int _basis_set_size;
            int _occupied_levels;
            int _unoccupied_levels;
            int _number_of_electrons;
            std::string _ECP;
            std::string _bsetype;


            std::map<int, std::vector<int> > _level_degeneracy;

            Eigen::VectorXd _mo_energies;
            Eigen::MatrixXd _mo_coefficients;

            Eigen::MatrixXd _overlap;
            Eigen::MatrixXd _vxc;

            std::vector< QMAtom* > _atoms;

            double _qm_energy;
            double _self_energy;

            Eigen::MatrixXd _mo_couplings;

            BasisSet _basis_set;

            // new variables for GW-BSE storage
            int _rpamin;
            int _rpamax;

            unsigned int _qpmin;
            unsigned int _qpmax;
            unsigned int _qptotal;

            unsigned int _bse_vmin;
            unsigned int _bse_vmax;
            unsigned int _bse_cmin;
            unsigned int _bse_cmax;
            unsigned int _bse_size;
            unsigned int _bse_vtotal;
            unsigned int _bse_ctotal;
            int _bse_nmax;

            double _ScaHFX;

            std::string _dftbasis;
            std::string _auxbasis;

            std::string _qm_package;

            // perturbative quasiparticle energies
            Eigen::MatrixXd _QPpert_energies;

            // quasiparticle energies and coefficients after diagonalization
            Eigen::VectorXd _QPdiag_energies;
            Eigen::MatrixXd _QPdiag_coefficients;
            // excitons



            MatrixXfd _eh_t;
            MatrixXfd _eh_s;
            VectorXfd _BSE_singlet_energies;
            MatrixXfd _BSE_singlet_coefficients;
            MatrixXfd _BSE_singlet_coefficients_AR;

            std::vector< tools::vec > _transition_dipoles;
            VectorXfd _BSE_triplet_energies;
            MatrixXfd _BSE_triplet_coefficients;

            Eigen::MatrixXd _BSE_singlet_couplings;
            Eigen::MatrixXd _BSE_triplet_couplings;
            int _couplingsA;
            int _couplingsB;



            std::vector< Eigen::VectorXd > _DqS_frag; // fragment charge changes in exciton

            std::vector< Eigen::VectorXd > _DqT_frag;

            Eigen::VectorXd _GSq_frag; // ground state effective fragment charges

            std::vector< Eigen::VectorXd > _popE_s;
            std::vector< Eigen::VectorXd > _popE_t;
            std::vector< Eigen::VectorXd > _popH_s;
            std::vector< Eigen::VectorXd > _popH_t;


            /**
             * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
             * @return A map with key as a level and a vector which is a list of close lying orbitals
             */
            bool CheckDegeneracy(double _energy_difference);

 
        };

    }
}


#endif /* __VOTCA_XTP_ORBITALS_H */
