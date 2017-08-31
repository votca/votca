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

#ifndef __VOTCA_XTP_ORBITALS_H
#define __VOTCA_XTP_ORBITALS_H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
//for openmp
#include <votca/xtp/votca_config.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/aobasis.h>
#include <votca/ctp/qmatom.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

#include <votca/ctp/logger.h>
#include <boost/format.hpp>
// Text archive that defines boost::archive::text_oarchive
// and boost::archive::text_iarchive
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// XML archive that defines boost::archive::xml_oarchive
// and boost::archive::xml_iarchive
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

// XML archive which uses wide characters (use for UTF-8 output ),
// defines boost::archive::xml_woarchive
// and boost::archive::xml_wiarchive
#include <boost/archive/xml_woarchive.hpp>
#include <boost/archive/xml_wiarchive.hpp>

// Binary archive that defines boost::archive::binary_oarchive
// and boost::archive::binary_iarchive
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include <votca/tools/constants.h>

#if (GWBSE_DOUBLE)
#define real_gwbse double
#else
#define real_gwbse float
#endif

namespace ub = boost::numeric::ublas;

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

            int getNumberOfElectrons() {
                return _number_of_electrons;
            };

            void setNumberOfElectrons(const int &electrons) {
                _number_of_electrons = electrons;
            }

            
            bool hasECP(){
                return ( _ECP !="") ? true : false;
            }
            
            string getECP() {
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
                return ( _overlap.size1() > 0) ? true : false;
            }

            const ub::symmetric_matrix<double> &AOOverlap() const {
                return _overlap;
            }

            ub::symmetric_matrix<double> &AOOverlap() {
                return _overlap;
            }

            // access to DFT molecular orbital energies, new, tested

            bool hasMOEnergies() {
                return ( _mo_energies.size() > 0) ? true : false;
            }

            const ub::vector<double> &MOEnergies() const {
                return _mo_energies;
            }

            ub::vector<double> &MOEnergies() {
                return _mo_energies;
            }

            // access to DFT molecular orbital energy of a specific level (in eV)

            double getEnergy(int level) {
                return ( hasMOEnergies()) ? votca::tools::conv::hrt2ev * _mo_energies[level - 1] : 0;
            }

            // access to DFT molecular orbital coefficients, new, tested

            bool hasMOCoefficients() {
                return ( _mo_coefficients.size1() > 0) ? true : false;
            }

            const ub::matrix<double> &MOCoefficients() const {
                return _mo_coefficients;
            }

            ub::matrix<double> &MOCoefficients() {
                return _mo_coefficients;
            }

            // access to DFT transfer integrals, new, tested

            bool hasMOCouplings() {
                return ( _mo_couplings.size1() > 0) ? true : false;
            }

            const ub::matrix<double> &MOCouplings() const {
                return _mo_couplings;
            }

            ub::matrix<double> &MOCouplings() {
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

            const std::vector< ctp::QMAtom* > &QMAtoms() const {
                return _atoms;
            }

            std::vector< ctp::QMAtom* > &QMAtoms() {
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
                return ( _vxc.size1() > 0) ? true : false;
            }

            ub::symmetric_matrix<double> &AOVxc() {
                return _vxc;
            }

            const ub::symmetric_matrix<double> &AOVxc() const {
                return _vxc;
            }

            // access to GW auxiliary basis set name

            bool hasGWbasis() {
                return ( !_gwbasis.empty()) ? true : false;
            }

            void setGWbasis(std::string basis) {
                _gwbasis = basis;
            }

            const std::string getGWbasis() const {
                return _gwbasis;
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

            void setBSEtype(string bsetype){_bsetype=bsetype;}
            string getBSEtype() const{return _bsetype;}
            
            
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
                for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
                    for (unsigned _c = 0; _c < _bse_ctotal; _c++) {
                        _index2v.push_back(_bse_vmin + _v);
                        _index2c.push_back(_bse_cmin + _c);
                    }
                }
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
                return ( _QPpert_energies.size1() > 0) ? true : false;
            }

            const ub::matrix<double> &QPpertEnergies() const {
                return _QPpert_energies;
            }

            ub::matrix<double> &QPpertEnergies() {
                return _QPpert_energies;
            }

            // access to diagonalized QP energies and wavefunctions

            bool hasQPdiag() {
                return ( _QPdiag_energies.size() > 0) ? true : false;
            }

            const ub::vector<double> &QPdiagEnergies() const {
                return _QPdiag_energies;
            }

            ub::vector<double> &QPdiagEnergies() {
                return _QPdiag_energies;
            }

            const ub::matrix<double> &QPdiagCoefficients() const {
                return _QPdiag_coefficients;
            }

            ub::matrix<double> &QPdiagCoefficients() {
                return _QPdiag_coefficients;
            }

            // access to eh interaction
            
         

            bool hasEHinteraction() {
                return ( _eh_d.size1() > 0) ? true : false;
            }

            const ub::matrix<real_gwbse> &eh_x() const {
                return _eh_x;
            }

            ub::matrix<real_gwbse> &eh_x() {
                return _eh_x;
            }

            const ub::matrix<real_gwbse> &eh_d() const {
                return _eh_d;
            }

            ub::matrix<real_gwbse> &eh_d() {
                return _eh_d;
            }

            // access to triplet energies and wave function coefficients

            bool hasBSETriplets() {
                return ( _BSE_triplet_energies.size() > 0) ? true : false;
            }

            const ub::vector<real_gwbse> &BSETripletEnergies() const {
                return _BSE_triplet_energies;
            }

            ub::vector<real_gwbse> &BSETripletEnergies() {
                return _BSE_triplet_energies;
            }

            const ub::matrix<real_gwbse> &BSETripletCoefficients() const {
                return _BSE_triplet_coefficients;
            }

            ub::matrix<real_gwbse> &BSETripletCoefficients() {
                return _BSE_triplet_coefficients;
            }

            // access to singlet energies and wave function coefficients

            bool hasBSESinglets() {
                return (_BSE_singlet_energies.size() > 0) ? true : false;
            }

            const ub::vector<real_gwbse> &BSESingletEnergies() const {
                return _BSE_singlet_energies;
            }

            ub::vector<real_gwbse> &BSESingletEnergies() {
                return _BSE_singlet_energies;
            }

            const ub::matrix<real_gwbse> &BSESingletCoefficients() const {
                return _BSE_singlet_coefficients;
            }

            ub::matrix<real_gwbse> &BSESingletCoefficients() {
                return _BSE_singlet_coefficients;
            }

            // for anti-resonant part in full BSE

            const ub::matrix<real_gwbse> &BSESingletCoefficientsAR() const {
                return _BSE_singlet_coefficients_AR;
            }

            ub::matrix<real_gwbse> &BSESingletCoefficientsAR() {
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
                return (_BSE_singlet_couplings.size1() > 0) ? true : false;
            }

            const ub::matrix<real_gwbse> &SingletCouplings() const {
                return _BSE_singlet_couplings;
            }

            ub::matrix<real_gwbse> &SingletCouplings() {
                return _BSE_singlet_couplings;
            }

            void setSingletCouplings(ub::matrix<real_gwbse> couplings) {
                _BSE_singlet_couplings = couplings;
            }

            // access to triplet coupling elements

            bool hasTripletCouplings() {
                return (_BSE_triplet_couplings.size1() > 0) ? true : false;
            }

            const ub::matrix<real_gwbse> &TripletCouplings() const {
                return _BSE_triplet_couplings;
            }

            ub::matrix<real_gwbse> &TripletCouplings() {
                return _BSE_triplet_couplings;
            }

            void setTripletCouplings(ub::matrix<real_gwbse> couplings) {
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
            ub::matrix<double> DensityMatrixGroundState();
            std::vector<ub::matrix<double> > DensityMatrixExcitedState(const string& spin,int state = 0);
            ub::matrix<double > TransitionDensityMatrix(const string& spin,int state = 0);


            double GetTotalEnergy(string _spintype, int _opt_state);

            // functions for analyzing fragment charges via Mulliken populations
            ub::vector<double> LoewdinPopulation(const ub::matrix<double>& _densitymatrix, const ub::matrix<double>& _overlapmatrix, int _frag);

            // access to fragment charges of singlet excitations

            bool hasFragmentChargesSingEXC() {
                return (_DqS_frag.size() > 0) ? true : false;
            }

            const std::vector< ub::vector<double> > &FragmentChargesSingEXC() const {
                return _DqS_frag;
            }

            std::vector< ub::vector<double> > &FragmentChargesSingEXC() {
                return _DqS_frag;
            }

            // access to fragment charges of triplet excitations

            bool hasFragmentChargesTripEXC() {
                return (_DqT_frag.size() > 0) ? true : false;
            }

            const std::vector< ub::vector<double> > &FragmentChargesTripEXC() const {
                return _DqT_frag;
            }

            std::vector< ub::vector<double> > &FragmentChargesTripEXC() {
                return _DqT_frag;
            }

            // access to fragment charges in ground state

            const ub::vector<double> &FragmentChargesGS() const {
                return _GSq_frag;
            }

            ub::vector<double> &FragmentChargesGS() {
                return _GSq_frag;
            }



            ub::vector<double> FragmentNuclearCharges(int _frag);




            // returns indeces of a re-sorted in a descending order vector of energies
            std::vector<int> SortEnergies();

            /** Adds a QM atom to the atom list */
            ctp::QMAtom* AddAtom(std::string _type,
                    double _x, double _y, double _z,
                    double _charge = 0, bool _from_environment = false) {
                ctp::QMAtom* pAtom = new ctp::QMAtom(_type, _x, _y, _z, _charge, _from_environment);
                _atoms.push_back(pAtom);
                return pAtom;
            }

            ctp::QMAtom* AddAtom(std::string _type, tools::vec pos,
                    double _charge = 0, bool _from_environment = false) {
                ctp::QMAtom* pAtom = new ctp::QMAtom(_type, pos, _charge, _from_environment);
                _atoms.push_back(pAtom);
                return pAtom;
            }

            ctp::QMAtom* AddAtom(ctp::QMAtom atom) {
                ctp::QMAtom* pAtom = new ctp::QMAtom(atom);
                _atoms.push_back(pAtom);
                return pAtom;
            }

            void setStorage(bool _store_orbitals, bool _store_overlap, bool _store_integrals) {
                // _has_mo_coefficients = _store_orbitals;
                //hasOverlap() = _store_overlap;
                // _has_integrals = _store_integrals;
                ;
            }

            void WritePDB(FILE *out, std::string tag = "");

            // reduces number of virtual orbitals to factor*number_of_occupied_orbitals
            void Trim(int factor);

            // reduces number of virtual orbitals to [HOMO-degG:LUMO+degL]
            void Trim(int degH, int degL);

            /** Loads orbitals from a file
             * Returns true if successful and does not throw an exception.
             * If exception is required, please use the << overload.
             */
            bool Load(std::string file_name);

            bool Save(std::string file_name);

            void LoadFromXYZ(std::string filename);
           

        private:
            std::vector<ub::matrix<double> > DensityMatrixExcitedState_R(const string& spin,int state = 0);
            std::vector<ub::matrix<double> >DensityMatrixExcitedState_AR(const string& spin,int state = 0);

            int _basis_set_size;
            int _occupied_levels;
            int _unoccupied_levels;
            int _number_of_electrons;
            string _ECP;
            
            string _bsetype;
            

            std::map<int, std::vector<int> > _level_degeneracy;

            ub::vector<double> _mo_energies;
            ub::matrix<double> _mo_coefficients;

            ub::symmetric_matrix<double> _overlap;
            ub::symmetric_matrix<double> _vxc;

            std::vector< ctp::QMAtom* > _atoms;

            double _qm_energy;
            double _self_energy;

            ub::matrix<double> _mo_couplings;

            bool _has_basis_set;
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
            std::string _gwbasis;

            std::string _qm_package;

            // perturbative quasiparticle energies
            ub::matrix<double> _QPpert_energies;

            // quasiparticle energies and coefficients after diagonalization
            ub::vector<double> _QPdiag_energies;
            ub::matrix<double> _QPdiag_coefficients;
            // excitons
            std::vector<int> _index2v;
            std::vector<int> _index2c;


            ub::matrix<real_gwbse> _eh_d;
            ub::matrix<real_gwbse> _eh_x;
            ub::vector<real_gwbse> _BSE_singlet_energies;
            ub::matrix<real_gwbse> _BSE_singlet_coefficients;
            ub::matrix<real_gwbse> _BSE_singlet_coefficients_AR;

            std::vector< tools::vec > _transition_dipoles;
            ub::vector<real_gwbse> _BSE_triplet_energies;
            ub::matrix<real_gwbse> _BSE_triplet_coefficients;

            ub::matrix<real_gwbse> _BSE_singlet_couplings;
            ub::matrix<real_gwbse> _BSE_triplet_couplings;
            int _couplingsA;
            int _couplingsB;



            std::vector< ub::vector<double> > _DqS_frag; // fragment charge changes in exciton

            std::vector< ub::vector<double> > _DqT_frag;

            ub::vector<double> _GSq_frag; // ground state effective fragment charges


        private:

            /**
             * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
             * @return A map with key as a level and a vector which is a list of close lying orbitals
             */
            bool CheckDegeneracy(double _energy_difference);

            // Allow serialization to access non-public data members
            friend class boost::serialization::access;

            //Allow  object to access non-public data members
            friend class Gaussian;
            friend class Turbomole;
            friend class NWChem;
            friend class Orca;
            friend class GW;

            // serialization itself (template implementation stays in the header)

            template<typename Archive>
            void serialize(Archive& ar, const unsigned int version) {
                //check with which votca version orbitals object was created
                if (version > 1) {
                    string test = "float";
#if (GWBSE_DOUBLE)
                    test = "double";
#endif

                    if (Archive::is_loading::value) {
                        string floatordouble = "float";
                        ar & floatordouble;
                        if (test != floatordouble) {
                            throw std::runtime_error((boost::format("This votca is compiled with %1%. The orbitals file you want to read in is compiled with %2%") % test % floatordouble).str());
                        }
                    } else {
                        ar & test;
                    }
                }
                ar & _basis_set_size;
                ar & _occupied_levels;
                ar & _unoccupied_levels;
                ar & _number_of_electrons;
                ar & _level_degeneracy;
                ar & _mo_energies;
                ar & _mo_coefficients;


                // symmetric matrix does not serialize by default
                if (Archive::is_saving::value) {
                    unsigned int size = _overlap.size1();
                    ar & size;
                }

                // copy the values back if loading
                if (Archive::is_loading::value) {
                    unsigned int size;
                    ar & size;
                    _overlap.resize(size);
                    if(version==0){
                        cerr<<endl;
                        cerr<<"WARNING! .orb file is of version 0. The overlap matrix will have the wrong ordering"<<endl;
                    }
                }

                for (unsigned int i = 0; i < _overlap.size1(); ++i)
                    for (unsigned int j = 0; j <= i; ++j)
                        ar & _overlap(i, j);

                ar & _atoms;
               
                
                ar & _qm_energy;
                ar & _qm_package;
                ar & _self_energy;
                ar & _mo_couplings;

                // GW-BSE storage
                if (version > 0) {
                    
                    

                    ar & _dftbasis;
                    ar & _gwbasis;
                    ar & _rpamin;
                    ar & _rpamax;
                    ar & _qpmin;
                    ar & _qpmax;
                    ar & _bse_vmin;
                    ar & _bse_vmax;
                    ar & _bse_cmin;
                    ar & _bse_cmax;
                    ar & _bse_nmax;
                    ar & _index2c;
                    ar & _index2v;
                    ar & _ScaHFX;
                    
                    if (Archive::is_loading::value && version < 4) {
                        _bsetype="TDA";
                        _ECP="ecp";   
                    } else {
                        ar & _bsetype;
                        ar & _ECP;
                    }
                    

                    if (Archive::is_loading::value && version < 3) {
                        ub::matrix<real_gwbse> temp;
                        ar &temp;
                        _QPpert_energies.resize(temp.size1(), temp.size2());
                        _QPpert_energies = 0.5 * temp;
                    } else {
                        ar & _QPpert_energies;
                    }

                    if (Archive::is_loading::value && version == 1) {
                        std::vector<double> temp;
                        ar &temp;
                        _QPdiag_energies.resize(temp.size());
                        for (unsigned i = 0; i < temp.size(); i++) {
                            _QPdiag_energies(i) = temp[i];
                        }
                    } else if (Archive::is_loading::value && version == 2) {
                        ub::vector<real_gwbse> temp;
                        ar &temp;
                        _QPdiag_energies.resize(temp.size());
                        _QPdiag_energies = 0.5 * temp;

                    } else {
                        ar & _QPdiag_energies;
                    }

                    ar & _QPdiag_coefficients;



                    if (Archive::is_loading::value && version < 3) {
                        ub::matrix<real_gwbse> temp;
                        ar &temp;
                        _eh_d.resize(temp.size1(), temp.size2());
                        _eh_d = 0.5 * temp;
                    } else {
                        ar & _eh_d;
                    }
                    if (Archive::is_loading::value && version < 3) {
                        ub::matrix<real_gwbse> temp;
                        ar &temp;
                        _eh_x.resize(temp.size1(), temp.size2());
                        _eh_x = 0.5 * temp;
                    } else {
                        ar & _eh_x;
                    }



                    if (Archive::is_loading::value && version == 1) {
                        std::vector<real_gwbse> temp;
                        ar &temp;
                        _BSE_singlet_energies.resize(temp.size());
                        for (unsigned i = 0; i < temp.size(); i++) {
                            _BSE_singlet_energies(i) = temp[i];
                        }
                    } else if (Archive::is_loading::value && version == 2) {
                        ub::vector<real_gwbse> temp;
                        ar &temp;
                        _BSE_singlet_energies.resize(temp.size());
                        _BSE_singlet_energies = 0.5 * temp;

                    } else {
                        ar & _BSE_singlet_energies;
                    }


                    ar & _BSE_singlet_coefficients;
                    
                    if (Archive::is_loading::value && version < 4) {
                        _BSE_singlet_coefficients_AR=ub::matrix<double>(0,0);
                    } else {
                        ar & _BSE_singlet_coefficients_AR;
                    }
                    

                    if (Archive::is_loading::value && version == 1) {
                        std::vector< std::vector<double> > temp;
                        ar &temp;
                        for (unsigned _i = 0; _i < temp.size(); _i++) {
                            ub::vector< double > vector_temp(3);
                            vector_temp(0) = temp[_i][0];
                            vector_temp(1) = temp[_i][1];
                            vector_temp(2) = temp[_i][2];
                            _transition_dipoles.push_back(vector_temp);
                        }
                    } else if (Archive::is_loading::value && version == 2) {
                        std::vector< ub::vector<double> > temp;
                        ar &temp;
                        for (unsigned _i = 0; _i < temp.size(); _i++) {
                            tools::vec vector_temp = tools::vec(temp[_i]);
                            _transition_dipoles.push_back(vector_temp);
                        }
                    } else {
                        ar & _transition_dipoles;
                    }



                    if (Archive::is_loading::value && version == 1) {
                        std::vector<real_gwbse> temp;
                        ar &temp;
                        _BSE_triplet_energies.resize(temp.size());
                        for (unsigned i = 0; i < temp.size(); i++) {
                            _BSE_triplet_energies(i) = 0.5 * temp[i];
                        }
                    } else if (Archive::is_loading::value && version == 2) {
                        ub::vector<real_gwbse> temp;
                        ar &temp;
                        _BSE_triplet_energies.resize(temp.size());
                        _BSE_triplet_energies = 0.5 * temp;
                    } else {
                        ar & _BSE_triplet_energies;
                    }

                    ar & _BSE_triplet_coefficients;


                    if (Archive::is_loading::value && version < 3) {
                        ub::matrix<real_gwbse> temp;
                        ar &temp;
                        _BSE_singlet_couplings.resize(temp.size1(), temp.size2());
                        _BSE_singlet_couplings = 0.5 * temp;
                    } else {
                        ar & _BSE_singlet_couplings;
                    }


                    if (Archive::is_loading::value && version < 3) {
                        ub::matrix<real_gwbse> temp;
                        ar &temp;
                        _BSE_triplet_couplings.resize(temp.size1(), temp.size2());
                        _BSE_triplet_couplings = 0.5 * temp;
                    } else {
                        ar & _BSE_triplet_couplings;
                    }



                    ar & _couplingsA;
                    ar & _couplingsB;

                    // symmetric matrix does not serialize by default
                    if (Archive::is_saving::value) {
                        unsigned int size = _vxc.size1();
                        ar & size;
                    }

                    // copy the values back if loading
                    if (Archive::is_loading::value) {
                        unsigned int size;
                        ar & size;
                        _vxc.resize(size);
                    }

                    for (unsigned int i = 0; i < _vxc.size1(); ++i)
                        for (unsigned int j = 0; j <= i; ++j)
                            ar & _vxc(i, j);

            if (Archive::is_loading::value && version < 4) {
                BasisSet _dftbasisset;
                _dftbasisset.LoadBasisSet(_dftbasis);

                if(!hasQMAtoms()){
                    throw runtime_error("Orbitals object has no QMAtoms");
                }
            AOBasis _dftbasis;
            _dftbasis.AOBasisFill(&_dftbasisset, QMAtoms());
            if(this->hasAOOverlap()){
                 _dftbasis.ReorderMatrix(_overlap,_qm_package , "xtp");

            }
            if(this->hasAOVxc()){
                
               
                if(_qm_package=="gaussian"){
                ub::matrix<double> vxc_full=_vxc;
                
                ub::matrix<double> _carttrafo=_dftbasis.getTransformationCartToSpherical(_qm_package);
                ub::matrix<double> _temp = ub::prod(_carttrafo, vxc_full);
                _vxc = ub::prod(_temp, ub::trans(_carttrafo));
                }
                 _dftbasis.ReorderMatrix(_vxc,_qm_package , "xtp");
                
            }   
            if(this->hasMOCoefficients()){
                _dftbasis.ReorderMOs(_mo_coefficients,_qm_package , "xtp");
            }  
                }

                } // end version 1: GW-BSE storage
            }// end of serialization
        };

    }
}

BOOST_CLASS_VERSION(votca::xtp::Orbitals, 4)

#endif /* __VOTCA_XTP_ORBITALS_H */

