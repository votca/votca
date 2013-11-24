/*
 *            Copyright 2009-2012 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/mbgft.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>
// #include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // MBPT MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        void MBGFT::CleanUp() {

        }
        
        
        /* 
         
         */
        
        bool MBGFT::Evaluate( Orbitals* _orbitals) {


            string _dft_package = _orbitals->getQMpackage();
            LOG(logDEBUG, *_pLog) << TimeStamp() << " DFT data was created by " << _dft_package << flush;

            // get atoms from orbitals object
            std::vector<QMAtom*>* _atoms = _orbitals->getAtoms();
            
            // reorder DFT data, load DFT basis set
            BasisSet dftbs;
            // string dftbasis_name("ubecppol");

            AOBasis dftbasis;

            dftbs.LoadBasisSet(_dftbasis_name);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

            //dftbasis.AOBasisFill(&dftbs, segments);
            dftbasis.AOBasisFill(&dftbs, *_atoms);

            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Basis of size " << dftbasis._AOBasisSize << flush;

   
            // convert _rpamax if needed 
            _homo = _orbitals->getNumberOfElectrons() -1 ; // indexed from 0
            _rpamin = 0; // lowest index occ min(gwa%mmin, screening%nsum_low) ! always 1
            if (_ranges == "default") {
                _rpamax = _orbitals->getNumberOfLevels() -1 ; // total number of levels
            } else if (_ranges == "factor") {
                _rpamax = _rpamaxfactor * _orbitals->getNumberOfLevels() -1; // total number of levels
            }
           
            // convert _qpmin and _qpmax if needed
            if (_ranges == "default") {
                _qpmin = 0;                     // indexed from 0
                _qpmax = 2 * _homo + 1;         // indexed from 0
            } else if (_ranges == "factor") {
                _qpmin = _orbitals->getNumberOfElectrons() - int( _qpminfactor * _orbitals->getNumberOfElectrons()) -1 ;
                _qpmax = _orbitals->getNumberOfElectrons() + int( _qpmaxfactor * _orbitals->getNumberOfElectrons()) -1 ;
            }else if (_ranges == "explicit") {
                _qpmin -= 1;
                _qpmax -= 1;
            }

            
           // set BSE band range indices 
            
            // anything else would be stupid!
            _bse_vmax = _homo;
            _bse_cmin = _homo +1;   
            
            if (_ranges == "default") {
                _bse_vmin = 0;                   // indexed from 0
                _bse_cmax = 2 * _homo + 1;       // indexed from 0
            } else if (_ranges == "factor") {
                _bse_vmin = _orbitals->getNumberOfElectrons() - int( _bseminfactor * _orbitals->getNumberOfElectrons()) -1 ;
                _bse_cmax = _orbitals->getNumberOfElectrons() + int( _bsemaxfactor * _orbitals->getNumberOfElectrons()) -1 ;
            }else if (_ranges == "explicit") {
                _bse_vmin -= 1;
                _bse_cmax -= 1;
            }
            _bse_vtotal = _bse_vmax - _bse_vmin +1 ;
            _bse_ctotal = _bse_cmax - _bse_cmin +1 ;
            _bse_size   = _bse_vtotal * _bse_ctotal;
            
            // indexing info
            for ( int _v = 0; _v < _bse_vtotal; _v++ ){
                for ( int _c = 0; _c < _bse_ctotal ; _c++){

                    _index2v.push_back( _bse_vmin + _v );
                    _index2c.push_back( _bse_cmin + _c );
                    
                }
            }

            
            // some QP - BSE consistency checks are required
            if ( _bse_vmin < _qpmin ) _qpmin = _bse_vmin;
            if ( _bse_cmax < _qpmax ) _qpmax = _bse_cmax;
            _qptotal = _qpmax - _qpmin +1 ;
            if ( _bse_nmax > _bse_size || _bse_nmax < 0 ) _bse_nmax = _bse_size;

            
                    
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Set RPA level range [" << _rpamin +1 << ":" << _rpamax +1 << "]" << flush;
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Set QP  level range [" << _qpmin +1 << ":" << _qpmax +1 << "]" << flush;
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Set BSE level range occ[" << _bse_vmin +1 << ":" << _bse_vmax +1 << "]  virt[" << _bse_cmin +1 << ":" << _bse_cmax +1 << "]" << flush;
                        
            
            // process the DFT data
            // a) form the expectation value of the XC functional in MOs
            ub::matrix<double> _dft_orbitals = *_orbitals->getOrbitals();
            // we have to do some cartesian -> spherical transformation for Gaussian
            ub::matrix<double> _vxc_ao;
            if ( _dft_package == "gaussian" ){

                ub::matrix<double> vxc_cart = (*_orbitals->getVxc()); 
                ub::matrix<double> _carttrafo;
                dftbasis.getTransformationCartToSpherical(_dft_package , _carttrafo );
                ub::matrix<double> _temp = ub::prod( _carttrafo, vxc_cart  );
                _vxc_ao = ub::prod( _temp, ub::trans( _carttrafo) );

            } else {
                _vxc_ao = (*_orbitals->getVxc()); 
            }

            // now get expectation values but only for those in _qpmin:_qpmax range
            ub::matrix<double> _mos = ub::project( _dft_orbitals ,  ub::range( _qpmin , _qpmax +1 ), ub::range(0, dftbasis._AOBasisSize ) );
            ub::matrix<double> _temp  = ub::prod( _vxc_ao, ub::trans( _mos ) ) ;
            _vxc  = ub::prod( _mos  , _temp );
            _vxc = 2.0 * _vxc;
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated exchange-correlation expectation values " << flush;
            

            
            // b) reorder MO coefficients depending on the QM package used to obtain the DFT data
            if (_dft_package != "votca") {
                // get reordering vector _dft_package -> Votca 
                vector<int> neworder;
                dftbasis.getReorderVector(_dft_package, neworder);
                // and reorder rows of _orbitals->_mo_coefficients() accordingly
                AOBasis::ReorderMOs(_dft_orbitals, neworder);
                // NWChem inverted sign for xz d-orbital
                if (_dft_package == "nwchem") {
                    // get vector with multipliers, e.g. NWChem -> Votca (bloody sign for d_xz)
                    vector<int> multiplier;
                    dftbasis.getMultiplierVector(_dft_package, multiplier);
                    // and reorder rows of _orbitals->_mo_coefficients() accordingly
                    AOBasis::MultiplyMOs(_dft_orbitals, multiplier);
                }
            }
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Converted DFT orbital coefficient order " << flush;

            // setting up ao_overlap_matrix
            list<string> elements;
            BasisSet gwbs;
            //string gwbasis_name("gwdefault");

            AOBasis gwbasis;
           // bool PPM_symmetric = true; // only PPM supported


            gwbs.LoadBasisSet(_gwbasis_name);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded GW Basis Set " << _gwbasis_name << flush;

            gwbasis.AOBasisFill(&gwbs, *_atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled GW Basis of size " << gwbasis._AOBasisSize << flush;

            
            
            // Testing momentum AOMatrix
            AOMomentum _dft_momentum;
            _dft_momentum.Initialize(dftbasis._AOBasisSize);
            _dft_momentum.Fill(&dftbasis);
  
            // now transition dipole elements for free interlevel transitions
           /* _interlevel_dipoles.resize(3);
            for ( int _i_comp = 0; _i_comp < 3; _i_comp++){
                // 1st: multiply each _dft_momentum component with range of empty states considered in BSE
                ub::matrix<double> _levels = ub::project( _dft_orbitals ,  ub::range( _bse_cmin , _bse_cmax +1 ), ub::range(0, dftbasis._AOBasisSize ) );
                _temp    = ub::prod( _dft_momentum._aomatrix[_i_comp] , ub::trans( _levels ) ) ;
                // 2nd: multiply this with range of occupied states considered in BSE
                _levels = ub::project( _dft_orbitals ,  ub::range( _bse_vmin , _bse_vmax +1 ), ub::range(0, dftbasis._AOBasisSize ) );
                _interlevel_dipoles[ _i_comp ] = ub::prod( _levels, _temp);
            }
            _temp.resize(0,0);
            _dft_momentum.Cleanup();
            LOG(logDEBUG, *pLog) << TimeStamp() << " Calculated free interlevel transition dipole moments (momentum) " << flush;         
            */
            
            // Testing electric dipole AOMatrix
            AODipole _dft_dipole;
            _dft_dipole.Initialize(dftbasis._AOBasisSize);
            _dft_dipole.Fill(&dftbasis);
  
            // now transition dipole elements for free interlevel transitions
            _interlevel_dipoles.resize(3);
            for ( int _i_comp = 0; _i_comp < 3; _i_comp++){
                // 1st: multiply each _dft_momentum component with range of empty states considered in BSE
                ub::matrix<double> _levels = ub::project( _dft_orbitals ,  ub::range( _bse_cmin , _bse_cmax +1 ), ub::range(0, dftbasis._AOBasisSize ) );
                _temp    = ub::prod( _dft_dipole._aomatrix[_i_comp] , ub::trans( _levels ) ) ;
                // 2nd: multiply this with range of occupied states considered in BSE
                _levels = ub::project( _dft_orbitals ,  ub::range( _bse_vmin , _bse_vmax +1 ), ub::range(0, dftbasis._AOBasisSize ) );
                _interlevel_dipoles[ _i_comp ] = ub::prod( _levels, _temp);
            }
            _temp.resize(0,0);
            _dft_dipole.Cleanup();
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated free interlevel transition dipole moments " << flush;          
            
            // get overlap matrix as AOOverlap
            AOOverlap _gwoverlap;
            // initialize overlap matrix
            _gwoverlap.Initialize(gwbasis._AOBasisSize);
            // Fill overlap
            _gwoverlap.Fill(&gwbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled GW Overlap matrix of dimension: " << _gwoverlap._aomatrix.size1() << flush;
            // _aooverlap.Print( "S" );

            // printing some debug info
            // _gwcoulomb.PrintIndexToFunction( &aobasis );

            
           
            
            // check eigenvalues of overlap matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            // _eigenvalues.resize(_gwoverlap._aomatrix.size1());
            // _eigenvectors.resize(_gwoverlap._aomatrix.size1(), _gwoverlap._aomatrix.size1());
            linalg_eigenvalues(_gwoverlap._aomatrix, _eigenvalues, _eigenvectors);
            // cout << _eigenvalues << endl;
            sort(_eigenvalues.begin(), _eigenvalues.end());
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of GW Overlap matrix : " << _eigenvalues[0] << flush;



            // get Coulomb matrix as AOCoulomb
            AOCoulomb _gwcoulomb;
            // initialize Coulomb matrix
            _gwcoulomb.Initialize(gwbasis._AOBasisSize);
            // Fill Coulomb matrix
            _gwcoulomb.Fill(&gwbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled GW Coulomb matrix of dimension: " << _gwcoulomb._aomatrix.size1() << flush;
            //_gwcoulomb.Print( "COU_in " );


            // PPM is symmetric, so we need to get the sqrt of the Coulomb matrix
            AOOverlap _gwoverlap_inverse;               // will also be needed in PPM itself
            AOOverlap _gwoverlap_cholesky_inverse;      // will also be needed in PPM itself
            _gwoverlap_inverse.Initialize( gwbasis._AOBasisSize);
            _gwcoulomb.Symmetrize(  _gwoverlap , gwbasis, _gwoverlap_inverse , _gwoverlap_cholesky_inverse );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Prepared GW Coulomb matrix for symmetric PPM " << flush;

            // calculate 3-center integrals,  convoluted with DFT eigenvectors

            // --- prepare a vector (gwdacay) of matrices (orbitals, orbitals) as container => M_mn
            // prepare 3-center integral object
            TCMatrix _Mmn;
            _Mmn.Initialize(gwbasis._AOBasisSize, _rpamin, _qpmax, _rpamin, _rpamax);
            _Mmn.Fill(gwbasis, dftbasis, _dft_orbitals);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated Mmn_beta (3-center-overlap x orbitals)  " << flush;
            //_Mmn.Print( "Mmn " );
            
            
            
            // for use in RPA, make a copy of _Mmn with dimensions (1:HOMO)(gwabasissize,LUMO:nmax)
            TCMatrix _Mmn_RPA;
            _Mmn_RPA.Initialize(gwbasis._AOBasisSize, _rpamin, _homo , _homo +1 , _rpamax);
            RPA_prepare_threecenters( _Mmn_RPA, _Mmn, gwbasis, _gwoverlap, _gwoverlap_inverse );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Prepared Mmn_beta for RPA  " << flush;

            // TODO: now, we can get rid of _gwoverlap_inverse
            // make _Mmn_RPA symmetric for use in RPA
            _Mmn_RPA.Symmetrize( _gwcoulomb._aomatrix  );
            //symmetrize_threecenters(_Mmn_RPA ,  _gwcoulomb._aomatrix);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Symmetrize Mmn_beta for RPA  " << flush;
            // _Mmn_RPA.Print( "Mmn_RPA" );
            
            // make _Mmn symmetric for use in self-energy calculation
            _Mmn.Symmetrize( _gwcoulomb._aomatrix  );
            //symmetrize_threecenters(_Mmn ,  _gwcoulomb._aomatrix);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Symmetrize Mmn_beta for self-energy  " << flush;


            // fix the frequencies for PPM
            _screening_freq = ub::zero_matrix<double>(2,2); // two frequencies
            //first one
            _screening_freq(0,0) = 0.0; // real part
            _screening_freq(0,1) = 0.0; // imaginary part
            //second one
            _screening_freq(1,0) = 0.0; // real part
            _screening_freq(1,1) = 1.0; // imaginary part

            ub::vector<double> _dft_energies = 2.0*(*_orbitals->getEnergies()); // getEnergies -> Hartree, we want Ryd
            
            // one entry to epsilon for each frequency
            _epsilon.resize( _screening_freq.size1() );
            
            // for symmetric PPM, we can initialize _epsilon with the overlap matrix!
            for ( int _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
                _epsilon[ _i_freq ] = _gwoverlap._aomatrix ;
            }
            
            
            
            
            // _gwoverlap is not needed further
            _gwoverlap.~AOOverlap();
            
            // determine epsilon from RPA
            RPA_calculate_epsilon( _Mmn_RPA, _screening_freq, _shift, _dft_energies );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated epsilon via RPA  " << flush;
  
            // _Mmn_RPA is not needed any further
            // _Mmn_RPA.~TCMatrix();
            _Mmn_RPA.Cleanup();
            
            // construct PPM parameters
            PPM_construct_parameters( _gwoverlap_cholesky_inverse._aomatrix );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Constructed PPM parameters  " << flush;
            
            // prepare threecenters for Sigma
            sigma_prepare_threecenters(  _Mmn );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Prepared threecenters for sigma  " << flush;

                        
            // calculate exchange part of sigma
            sigma_x_setup( _Mmn );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated exchange part of Sigma  " << flush;
            
            // TOCHECK: get rid of _ppm_phi?
  
            
            // calculate correlation part of sigma
            sigma_c_setup( _Mmn, _dft_energies   );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated correlation part of Sigma  " << flush;
 
            /* One could also save the qp_energies directly to the orbitals object...
                       // now copy energies in map into orbitals object matrix
            (_orbitals->_QPpert_energies).resize( _levels, 5 );
            for (size_t itl=0; itl < _levels; itl++){
               for (size_t ite=0; ite<5; ite++) {
                   _orbitals->_QPpert_energies(itl,ite) = _energies[itl+1][ite];
               }
            } 
             */
            
            
            
            // Output of quasiparticle energies after all is done:
            _pLog->setPreface(logINFO, "\n");
            
            LOG(logINFO,*_pLog) << (format("  ====== Perturbative quasiparticle energies (Rydberg) ====== ")).str() << flush;
                        
            for ( int _i = 0 ; _i < _qptotal ; _i++ ){
                if ( (_i + _qpmin) == _homo ){
                    LOG(logINFO,*_pLog) << (format("  HOMO  = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = %4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") % (_i+_qpmin+1) % _dft_energies( _i + _qpmin ) % _vxc(_i,_i) % _sigma_x(_i,_i) % _sigma_c(_i,_i) % _qp_energies(_i) ).str() << flush;
                } else if ( (_i + _qpmin) == _homo+1 ){
                    LOG(logINFO,*_pLog) << (format("  LUMO  = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = %4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") % (_i+_qpmin+1) % _dft_energies( _i + _qpmin ) % _vxc(_i,_i) % _sigma_x(_i,_i) % _sigma_c(_i,_i) % _qp_energies(_i) ).str() << flush;                    
                    
                }else {
                LOG(logINFO,*_pLog) << (format("  Level = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = %4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") % (_i+_qpmin+1) % _dft_energies( _i + _qpmin ) % _vxc(_i,_i) % _sigma_x(_i,_i) % _sigma_c(_i,_i) % _qp_energies(_i) ).str() << flush;
                }
            }

            // constructing full quasiparticle Hamiltonian and diagonalize, if requested
            if ( _do_qp_diag || _do_bse_singlets || _do_bse_triplets ){
            FullQPHamiltonian();
            if ( _do_qp_diag ) {
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Full quasiparticle Hamiltonian  " << flush;
            LOG(logINFO, *_pLog) << (format("  ====== Diagonalized quasiparticle energies (Rydberg) ====== ")).str() << flush;
            for (int _i = 0; _i < _qptotal; _i++) {
                if ((_i + _qpmin) == _homo) {
                    LOG(logINFO, *_pLog) << (format("  HOMO  = %1$4d PQP = %2$+1.4f DQP = %3$+1.4f ") % (_i + _qpmin + 1) % _qp_energies(_i) % _qp_diag_energies(_i)).str() << flush;
                } else if ((_i + _qpmin) == _homo + 1) {
                    LOG(logINFO, *_pLog) << (format("  LUMO  = %1$4d PQP = %2$+1.4f DQP = %3$+1.4f ") % (_i + _qpmin + 1) % _qp_energies(_i) % _qp_diag_energies(_i)).str() << flush;

                } else {
                    LOG(logINFO, *_pLog) << (format("  Level = %1$4d PQP = %2$+1.4f DQP = %3$+1.4f ") % (_i + _qpmin + 1) % _qp_energies(_i) % _qp_diag_energies(_i)).str() << flush;
                }
            }
            }
            }
            // proceed only if BSE requested
            if ( _do_bse_singlets || _do_bse_triplets ){
            BSE_qp_setup();
            // constructing electron-hole interaction for BSE
            if ( _do_bse_singlets ){
            // calculate exchange part of eh interaction, only needed for singlets
            BSE_x_setup(  _Mmn );
            LOG(logINFO,*_pLog) << TimeStamp() << " Exchange part of e-h interaction " << flush; 
            }
            
            
            // calculate direct part of eh interaction, needed for singlets and triplets
            BSE_d_setup ( _Mmn );
            LOG(logINFO,*_pLog) << TimeStamp() << " Direct part of e-h interaction " << flush; 


            if ( _do_bse_triplets ){
            BSE_solve_triplets();
            LOG(logINFO,*_pLog) << TimeStamp() << " Solved BSE for triplets " << flush; 
            
            // expectation values, contributions from e-h coupling
            std::vector<double> _contrib_x(10,0.0);
            std::vector<double> _contrib_d(10,0.0);
            std::vector<double> _contrib_qp(10,0.0);
            for (int _i_exc = 0; _i_exc < 10; _i_exc++){
                
                ub::matrix<double> _slice = ub::project( _bse_triplet_coefficients ,  ub::range( 0 , _bse_size ), ub::range(_i_exc, _i_exc + 1 ) );
                
                ub::matrix<double> _temp = ub::prod( _eh_d, _slice );
                ub::matrix<double>  _res =  ub::prod( ub::trans( _slice) , _temp  );
                _contrib_d[_i_exc] = -_res(0,0);
                
                _temp = ub::prod( _eh_qp, _slice );
                _res =  ub::prod( ub::trans( _slice) , _temp  );
                _contrib_qp[_i_exc] = _res(0,0);
                
            }
            
            
            
            LOG(logINFO, *_pLog) << (format("  ====== 10 lowest triplet energies (eV) ====== ")).str() << flush;
            for (int _i = 0; _i < 10; _i++) {
                LOG(logINFO, *_pLog) << (format("  T = %1$4d Omega = %2$+1.4f <FT> = %3$+1.4f <K_x> = %4$+1.4f <K_d> = %5$+1.4f") % (_i + 1) % (13.6058 * _bse_triplet_energies( _i )) % (13.6058 * _contrib_qp[_i]) % (13.6058 * _contrib_x[_i]) % (13.6058 * _contrib_d[ _i ]) ).str() << flush;

                for ( int _i_bse = 0 ; _i_bse < _bse_size; _i_bse++){
                    // if contribution is larger than 0.2, print
                    double _weight = pow(_bse_triplet_coefficients( _i_bse, _i),2);
                    if ( _weight > 0.2 ){
                        LOG(logINFO, *_pLog) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") % (_homo-_index2v[_i_bse]) % ( _index2c[_i_bse] - _homo -1) % (100.0 * _weight) ).str() << flush;
                    }
                }
                LOG(logINFO, *_pLog) << (format("   ")).str() << flush;
            }
            } // do_triplets
            
            if ( _do_bse_singlets ){
            BSE_solve_singlets();
            LOG(logINFO,*_pLog) << TimeStamp() << " Solved BSE for singlets " << flush; 


            // expectation values, contributions from e-h coupling
                             std::vector<double> _contrib_x(10,0.0);
            std::vector<double> _contrib_d(10,0.0);
            std::vector<double> _contrib_qp(10,0.0);           
            for (int _i_exc = 0; _i_exc < 10; _i_exc++){

                ub::matrix<double> _slice = ub::project( _bse_singlet_coefficients ,  ub::range( 0 , _bse_size ), ub::range(_i_exc, _i_exc + 1 ) );
                
                ub::matrix<double> _temp = ub::prod( _eh_x, _slice );
                ub::matrix<double>  _res =  ub::prod( ub::trans( _slice) , _temp  );
                _contrib_x[_i_exc] = 2.0 * _res(0,0);

                _temp = ub::prod( _eh_d, _slice );
                _res =  ub::prod( ub::trans( _slice) , _temp  );
                _contrib_d[_i_exc] = -_res(0,0);
                
                _temp = ub::prod( _eh_qp, _slice );
                _res =  ub::prod( ub::trans( _slice) , _temp  );
                _contrib_qp[_i_exc] = _res(0,0);
                
            }
            
            // for transition dipole moments
            std::vector<std::vector<double> > _transition_dipoles;
            std::vector<double> _oscillator_strength;
            std::vector<double> _transition_dipole_strength;
            for (int _i_exc = 0; _i_exc < 10 ; _i_exc++){
                std::vector<double> _tdipole(3,0.0);
                
                for ( int _v = 0 ; _v < _bse_vtotal ; _v++ ){
                    for ( int _c = 0; _c < _bse_ctotal ; _c++){
 
                        int index_vc = _bse_ctotal * _v + _c;
                    
                        for ( int _i_comp =0 ; _i_comp < 3; _i_comp++){
                                _tdipole[ _i_comp ] += _bse_singlet_coefficients( index_vc , _i_exc   ) * _interlevel_dipoles[_i_comp](_v,_c);
                        }
                        
                    }
                    
                }
		//                cout << _i_exc << " : " << _tdipole[0] << " : "  << _tdipole[1] << " : " << _tdipole[2] << endl;
                _transition_dipoles.push_back( _tdipole );
                _transition_dipole_strength.push_back( (_tdipole[0]*_tdipole[0] + _tdipole[1]*_tdipole[1] + _tdipole[2]*_tdipole[2]  )  );               
                _oscillator_strength.push_back( _transition_dipole_strength[_i_exc] * 1.0/3.0 * (_bse_singlet_energies( _i_exc )) ) ;
            }
            
            
            LOG(logINFO, *_pLog) << (format("  ====== 10 lowest singlet energies (eV) ====== ")).str() << flush;
            for (int _i = 0; _i < 10; _i++) {

                LOG(logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.4f <FT> = %3$+1.4f <K_x> = %4$+1.4f <K_d> = %5$+1.4f") % (_i + 1) % (13.6058 * _bse_singlet_energies( _i )) % (13.6058 * _contrib_qp[_i]) % (13.6058 * _contrib_x[_i]) % (13.6058 * _contrib_d[ _i ]) ).str() << flush;
                LOG(logINFO, *_pLog) << (format("           TrDipole length gauge   dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") %  (_transition_dipoles[_i][0]) % (_transition_dipoles[_i][1]) % (_transition_dipoles[_i][2]) % (_transition_dipole_strength[_i]) %  (_oscillator_strength[_i])).str() << flush;
                for ( int _i_bse = 0 ; _i_bse < _bse_size; _i_bse++){
                    // if contribution is larger than 0.2, print
                    double _weight = pow(_bse_triplet_coefficients( _i_bse, _i),2);
                    if ( _weight > 0.2 ){
                        LOG(logINFO, *_pLog) << (format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") % (_homo-_index2v[_i_bse]) % ( _index2c[_i_bse] - _homo -1) % (100.0 * _weight) ).str() << flush;
                    }
                }
                LOG(logINFO, *_pLog) << (format("   ")).str() << flush;
            }
            
            }
            }
            

            
            return true;
        }


        
        void MBGFT::BSE_qp_setup(){
            _eh_qp = ub::zero_matrix<double>( _bse_size , _bse_size );
            
            
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _eh_qp( _index_vc , _index_vc ) += _vxc(_c1 + _bse_cmin ,_c1 + _bse_cmin ) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _eh_qp( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_cmin ,_c2 + _bse_cmin );
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _eh_qp( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            
        }
        
        
        
        void MBGFT::BSE_solve_triplets(){
            
            ub::matrix<double> _bse = _eh_qp -_eh_d;
            /*
            // add full QP Hamiltonian contributions to free transitions
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _bse( _index_vc , _index_vc ) += _vxc(_c1 + _bse_cmin ,_c1 + _bse_cmin ) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _bse( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_cmin ,_c2 + _bse_cmin );
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _bse( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            */
            
            linalg_eigenvalues( _bse, _bse_triplet_energies, _bse_triplet_coefficients, _bse_nmax);
        }
        
        void MBGFT::BSE_solve_singlets(){
            
            ub::matrix<double> _bse = -_eh_d + 2.0 * _eh_x;


            // add full QP Hamiltonian contributions to free transitions
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _bse( _index_vc , _index_vc ) += _vxc(_c1 + _bse_cmin ,_c1 + _bse_cmin) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _bse( _index_vc , _index_vc2 ) += _vxc(_c1 + _bse_cmin ,_c2 + _bse_cmin);
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _bse( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            // _bse_singlet_energies.resize(_bse_singlet_coefficients.size1());
            linalg_eigenvalues(_bse, _bse_singlet_energies, _bse_singlet_coefficients, _bse_nmax);
            
          //  cout << TimeStamp() << " with GSL " << endl;
        }
        
        
        void MBGFT::BSE_d_setup (const TCMatrix& _Mmn){
            // gwbasis size
            size_t _gwsize = _Mmn[0].size1();

            // messy procedure, first get two matrices for occ and empty subbparts
            // store occs directly transposed
            ub::matrix<double> _storage_v = ub::zero_matrix<double>(  _bse_vtotal * _bse_vtotal , _gwsize );
            #pragma omp parallel for
            for ( size_t _v1 = 0; _v1 < _bse_vtotal; _v1++){
                const ub::matrix<double>& Mmn = _Mmn[_v1];
                for ( size_t _v2 = 0; _v2 < _bse_vtotal; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_v( _index_vv , _i_gw ) = Mmn( _i_gw , _v2  );
                    }
                }
            }
            
            
            ub::matrix<double> _storage_c = ub::zero_matrix<double>( _gwsize, _bse_ctotal * _bse_ctotal );
            #pragma omp parallel for
            for ( size_t _c1 = 0; _c1 < _bse_ctotal; _c1++){
                const ub::matrix<double>& Mmn = _Mmn[_c1 + _bse_cmin];
                for ( size_t _c2 = 0; _c2 < _bse_ctotal; _c2++){
                    size_t _index_cc = _bse_ctotal * _c1 + _c2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_c( _i_gw , _index_cc ) = Mmn( _i_gw , _c2 + _bse_cmin );
                    }
                }
            }
            
            // store elements in a vtotal^2 x ctotal^2 matrix
            // cout << "BSE_d_setup 1 [" << _storage_v.size1() << "x" << _storage_v.size2() << "]\n" << std::flush;
            ub::matrix<double> _storage_prod = ub::prod( _storage_v , _storage_c );
            
            
            
            // now patch up _storage for screened interaction
            #pragma omp parallel for
            for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){  
                double _ppm_factor = sqrt( _ppm_weight( _i_gw ));
                for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_vtotal) ; _v++){
                    _storage_v( _v , _i_gw ) = _ppm_factor * _storage_v(_v , _i_gw );
                }
                for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_ctotal) ; _c++){
                    _storage_c( _i_gw , _c ) = _ppm_factor * _storage_c( _i_gw , _c  );
                }
            }
            
            // multiply and subtract from _storage_prod
            // cout << "BSE_d_setup 2 [" << _storage_c.size1() << "x" << _storage_c.size2() << "]\n" << std::flush;
            _storage_prod -= ub::prod( _storage_v , _storage_c );
            
            
            // finally resort into _eh_d and multiply by to for Rydbergs
            // can be limited to upper diagonal !
            _eh_d = ub::zero_matrix<double>( _bse_size , _bse_size );
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_vc1 = _bse_ctotal * _v1 + _c1 ;
                              
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_vc2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_cc  = _bse_ctotal * _c1 + _c2;

                            _eh_d( _index_vc1 , _index_vc2 ) = 2.0 * _storage_prod( _index_vv , _index_cc ); 

                            
                        }
                    }
                }
            }
            
            
        }
        
        
        
        void MBGFT::BSE_x_setup(const TCMatrix& _Mmn){
            
            /* unlike the fortran code, we store eh interaction directly in
             * a suitable matrix form instead of a four-index array
             */
            
            // allocate eh_x
            _eh_x = ub::zero_matrix<double>( _bse_size , _bse_size );
            
            // gwbasis size
            size_t _gwsize = _Mmn[0].size1();
            
            // get a different storage for 3-center integrals we need
            ub::matrix<double> _storage = ub::zero_matrix<double>( _gwsize , _bse_size);
            
            // occupied levels
            #pragma omp parallel for
            for ( size_t _v = 0; _v < _bse_vtotal ; _v++ ){
                const ub::matrix<double>& Mmn = _Mmn[_v];
                // empty levels
                for (size_t _c =0 ; _c < _bse_ctotal ; _c++ ){
                    size_t _index_vc = _bse_ctotal * _v + _c ;
                    for (size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){
                        _storage( _i_gw, _index_vc ) = Mmn( _i_gw, _c + _bse_cmin);
                    }
                }
            }
            
            // with this storage, _eh_x is obtained by matrix multiplication
            _eh_x = ub::prod( ub::trans( _storage ), _storage );
            _eh_x = 2.0 * _eh_x; // Rydberg
  
            
        }
        
        
        

        void MBGFT::FullQPHamiltonian(){
            
            // constructing full QP Hamiltonian, storage in vxc
            _vxc = -_vxc + _sigma_x + _sigma_c;
            // diagonal elements are given by _qp_energies
            for (int _m = 0; _m < _vxc.size1(); _m++ ){
              _vxc( _m,_m ) = _qp_energies( _m );
            }
                   
            
            
            if ( _do_qp_diag ){
            /* diagonalize, since _vxc will be needed in BSE, and GSL
             * destroys the input array, we need to make a local copy first
             */
            
            // get eigenvalues and eigenvectors of this matrix
            ub::matrix<double> _temp = _vxc;
            _qp_diag_energies.resize(_temp.size1());
            _qp_diag_coefficients.resize(_temp.size1(), _temp.size1());
            linalg_eigenvalues(_temp, _qp_diag_energies, _qp_diag_coefficients);

            // TODO storage -> orbitals

            }
            
        }
        
        
        void MBGFT::sigma_c_setup(const TCMatrix& _Mmn, const ub::vector<double>& _edft){
            
            // iterative refinement of qp energies
            int _max_iter = 5;
            int _levelsum = _Mmn[0].size2(); // total number of bands
            int _gwsize  = _Mmn[0].size1(); // size of the GW basis
            const double pi = boost::math::constants::pi<double>();
            

            // initial _qp_energies are dft energies
            _qp_energies = _edft; // RANGES!
            
     

	    // only diagonal elements except for in final iteration
            for ( int _i_iter = 0 ; _i_iter < _max_iter-1 ; _i_iter++ ){
                
	      // initialize sigma_c to zero at the beginning of each iteration
	      _sigma_c = ub::zero_matrix<double>(_qptotal,_qptotal);

	      // loop over all GW levels
              #pragma omp parallel for
	      for (int _gw_level = 0; _gw_level < _qptotal ; _gw_level++ ){
                  
                const ub::matrix<double>& Mmn = _Mmn[ _gw_level ];
              
		// loop over all functions in GW basis
		for ( int _i_gw = 0; _i_gw < _gwsize ; _i_gw++ ){
                    
		  // loop over all bands
		  for ( int _i = 0; _i < _levelsum ; _i++ ){
                    
		    double occ = 1.0;
		    if ( _i > _homo ) occ = -1.0; // sign for empty levels
                                                    
		    // energy denominator
		    double _denom = _qp_energies( _gw_level ) - _qp_energies( _i ) + occ*_ppm_freq( _i_gw );
                    
		    double _stab = 1.0;
		    if ( std::abs(_denom) < 0.5 ) {
		      _stab = 0.5 * ( 1.0 - cos(2.0 * pi * std::abs(_denom) ) );
		    }
                            
		    double _factor = _ppm_weight( _i_gw ) * _ppm_freq( _i_gw) * _stab/_denom; // contains conversion factor 2!
		    
		    // sigma_c diagonal elements
		    _sigma_c( _gw_level , _gw_level ) += _factor * Mmn(_i_gw, _i) *  Mmn(_i_gw, _i);  
                            
		  }// bands
                        
		}// GW functions
		
		// update _qp_energies
		_qp_energies( _gw_level ) = _edft( _gw_level ) + _sigma_x(_gw_level, _gw_level) + _sigma_c(_gw_level,_gw_level) - _vxc(_gw_level,_gw_level);
                    
	      }// all bands
                
            } // iterations
            

	    // in final step, also calc offdiagonal elements
	    // initialize sigma_c to zero at the beginning
	    _sigma_c = ub::zero_matrix<double>(_qptotal,_qptotal);

            // loop over col  GW levels
            #pragma omp parallel for
	    for (int _gw_level = 0; _gw_level < _qptotal ; _gw_level++ ){
              
                const ub::matrix<double>& Mmn =  _Mmn[ _gw_level ];

		// loop over all functions in GW basis
		for ( int _i_gw = 0; _i_gw < _gwsize ; _i_gw++ ){
                    

		  // loop over all screening levels
		  for ( int _i = 0; _i < _levelsum ; _i++ ){
                    
		    double occ = 1.0;
		    if ( _i > _homo ) occ = -1.0; // sign for empty levels
                    
		    // energy denominator
		    double _denom = _qp_energies( _gw_level ) - _qp_energies( _i ) + occ*_ppm_freq( _i_gw );
                    
		    double _stab = 1.0;
		    if ( std::abs(_denom) < 0.5 ) {
		      _stab = 0.5 * ( 1.0 - cos(2.0 * pi * std::abs(_denom) ) );
		    }
                    
		    double _factor = _ppm_weight( _i_gw ) * _ppm_freq( _i_gw) *   Mmn(_i_gw, _i) *_stab/_denom; // contains conversion factor 2!
		    
		    // loop over row GW levels
		    for ( int _m = 0 ; _m < _qptotal ; _m++) {

		    
		    // sigma_c all elements
		    _sigma_c( _m , _gw_level ) += _factor * _Mmn[_m](_i_gw, _i) ;  //_submat(_i_gw,_i);
	                      
		  }// screening levels
		}// GW functions
	      }// GW row 
	      _qp_energies( _gw_level ) = _edft( _gw_level ) + _sigma_x(_gw_level,_gw_level) + _sigma_c(_gw_level,_gw_level) - _vxc(_gw_level,_gw_level);
	    } // GW col
    
        } // sigma_c_setup

        void MBGFT::sigma_x_setup(const TCMatrix& _Mmn){
        
            // initialize sigma_x
            _sigma_x = ub::zero_matrix<double>(_qptotal,_qptotal);
            int _size  = _Mmn[0].size1();

            // band 1 loop over all GW levels
            #pragma omp parallel for
            for ( int _m1 = 0 ; _m1 < _qptotal ; _m1++ ){
                
                const ub::matrix<double>& M1mn =  _Mmn[ _m1 ];
                
                // band 2 loop over all GW levels
                for ( int _m2 = 0 ; _m2 < _qptotal ; _m2++ ){
                    
                    const ub::matrix<double>& M2mn =  _Mmn[ _m2 ];
                    
                    // loop over all basis functions
                    for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                        // loop over all occupied bands used in screening
                        for ( int _i_occ = 0 ; _i_occ <= _homo ; _i_occ++ ){
                            _sigma_x( _m1, _m2 ) -= 2.0 * M1mn( _i_gw , _i_occ ) * M2mn( _i_gw , _i_occ );
                        } // occupied bands
                    } // gwbasis functions
                } // level 2
            } // level 1
        
        }



        

        void MBGFT::sigma_prepare_threecenters(TCMatrix& _Mmn){
            #pragma omp parallel for
            for ( int _m_level = 0 ; _m_level < _Mmn.get_mtot(); _m_level++ ){
                // get Mmn for this _m_level
                // ub::matrix<double> _temp = ub::trans(  _Mmn.matrix()( _m_level )   );
                // and multiply with _ppm_phi = eigenvectors of epsilon
                // POTENTIAL BUG
                _Mmn[ _m_level ] = ub::prod(  _ppm_phi , _Mmn[ _m_level ] );
            }
        }        
        
        void MBGFT::PPM_construct_parameters(  ub::matrix<double>& _overlap_cholesky_inverse ){
            
            // multiply with L-1^t from the right
            ub::matrix<double> _overlap_cholesky_inverse_transposed = ub::trans( _overlap_cholesky_inverse );
            ub::matrix<double> _temp = ub::prod( _epsilon[0] , _overlap_cholesky_inverse_transposed );
            // multiply with L-1 from the left
            _temp = ub::prod( _overlap_cholesky_inverse, _temp );
            
            // get eigenvalues and eigenvectors of this matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            _eigenvalues.resize(_temp.size1());
            _eigenvectors.resize(_temp.size1(), _temp.size1());
            linalg_eigenvalues(_temp, _eigenvalues, _eigenvectors);
            
            // multiply eigenvectors with overlap_cholesky_inverse_transpose and store as eigenvalues of epsilon
            _ppm_phi = ub::prod( _overlap_cholesky_inverse_transposed , _eigenvectors ); 

 
            
            // store PPM weights from eigenvalues
            _ppm_weight.resize( _eigenvalues.size() );
            for ( int _i = 0 ; _i <  _eigenvalues.size(); _i++   ){
                _ppm_weight(_i) = 1.0 - 1.0/_eigenvalues(_i);
            }

            // determine PPM frequencies
            _ppm_freq.resize( _eigenvalues.size() );
            // a) phi^t * epsilon(1) * phi 
            _temp = ub::prod( ub::trans( _ppm_phi ) , _epsilon[1] );
            _eigenvectors  = ub::prod( _temp ,  _ppm_phi  );
            // b) invert
            _temp = ub::zero_matrix<double>( _eigenvalues.size(),_eigenvalues.size() )  ;
            linalg_invert( _eigenvectors , _temp ); //eigenvectors is destroyed after!
            // c) PPM parameters -> diagonal elements
            for ( int _i = 0 ; _i <  _eigenvalues.size(); _i++   ){
                
                double _nom  = _temp( _i, _i ) - 1.0;
                
                // only purely imaginary frequency assumed
                if ( _screening_freq(1,0) != 0.0 ) {
                    cerr << " mixed frequency! real part: " << _screening_freq( 1, 0 ) << " imaginary part: "  << _screening_freq( 1 , 1 ) << flush;
                    exit(1);
                } else {
                    
                    double _frac = -1.0 * _nom/(_nom + _ppm_weight( _i )) * _screening_freq(1,1) * _screening_freq(1,1) ;
                    _ppm_freq ( _i ) =  sqrt( std::abs(_frac )) ;

		    if ( _ppm_weight(_i) < 1.e-5 ){
		      _ppm_weight(_i) = 0.0;
		      _ppm_freq(_i)   = 1.0;

		    }

                }

            }
            
            // will be needed transposed later
            _ppm_phi = ub::trans( _ppm_phi );
            
            // epsilon can be deleted
            _epsilon[0].resize(0,0);
            _epsilon[1].resize(0,0);
                   
            
        }
        
        

        void MBGFT::RPA_calculate_epsilon(TCMatrix& _Mmn_RPA, ub::matrix<double> _screening_freq, double _shift, ub::vector<double>& _dft_energies){
            
            int _size = _Mmn_RPA[0].size1(); // size of gwbasis
            
            // loop over frequencies
            for ( int _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
                
                // loop over occupied levels -> vector index of _Mmn_RPA
                // by default all variable shared, except for one defined in parallel region
                #pragma omp parallel for 
                for ( int _m_level = 0; _m_level < _Mmn_RPA.get_mtot() ; _m_level++ ){
                    int index_m = _Mmn_RPA.get_mmin();
                    const ub::matrix<double>& Mmn_RPA =  _Mmn_RPA[ _m_level ];

                    // a temporary matrix, that will get filled in empty levels loop
                    ub::matrix<double> _temp = ub::zero_matrix<double>( _Mmn_RPA.get_ntot(), _size );
                    
                        
                    // loop over empty levels
                    for ( int _n_level = 0 ; _n_level < _Mmn_RPA.get_ntot() ; _n_level++ ){
                        int index_n = _Mmn_RPA.get_nmin();
                        
                        
                        double _deltaE = _shift + _dft_energies( _n_level + index_n ) - _dft_energies( _m_level + index_m ); // get indices and units right!!!
                        double _energy_factor;
                        // this only works, if we have either purely real or purely imaginary frequencies
                        if ( _screening_freq( _i_freq, 0) == 0.0 ) {
                            // purely imaginary
                            _energy_factor = 8.0 * _deltaE / (_deltaE*_deltaE + _screening_freq( _i_freq, 1) *  _screening_freq( _i_freq, 1 ));
                        } else if ( _screening_freq( _i_freq, 1) == 0.0  ) {
                            // purely real
                            _energy_factor = 4.0 * (1.0 / (_deltaE - _screening_freq( _i_freq, 0 ) ) +  1.0 / (_deltaE + _screening_freq( _i_freq, 0 ) ) );
                        } else {
                            // mixed -> FAIL
                            cerr << " mixed frequency! real part: " << _screening_freq( _i_freq, 0 ) << " imaginary part: "  << _screening_freq( _i_freq, 1 ) << flush;
                            exit(1);
                        }

                        
                        
                        // _temp = _energy_factor * ub::trans( Mmn_RPA );
                        for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                            _temp( _n_level , _i_gw ) = _energy_factor * Mmn_RPA( _i_gw , _n_level );
                         } // matrix size
                        
                    } // empty levels

                   // now multiply and add to epsilon
                    ub::matrix<double> _add = ub::prod( Mmn_RPA , _temp  );
                   #pragma omp critical
                    {
                   _epsilon[ _i_freq ] += _add;// ub::prod( Mmn_RPA , _temp  );
                    }
                } // occupied levels
                
            } // loop over frequencies
            
            
        }
        
        
   
    
    void MBGFT::RPA_prepare_threecenters( TCMatrix& _Mmn_RPA, TCMatrix& _Mmn_full, AOBasis& gwbasis, AOMatrix& gwoverlap, AOMatrix& gwoverlap_inverse     ){
        
         
        // loop over m-levels in _Mmn_RPA
        #pragma omp parallel for
        for ( int _m_level = 0; _m_level < _Mmn_RPA.size() ; _m_level++ ){
        
            ub::matrix<double> _temp = ub::prod( gwoverlap_inverse._aomatrix , _Mmn_full[ _m_level ] );

            // loop over n-levels in _Mmn_full 
            for ( int _n_level = 0; _n_level < _Mmn_full.get_ntot() ; _n_level++ ){

                double sc_plus  = 0.0;
                double sc_minus = 0.0;
                
                // loop over gwbasis shells
                for (vector< AOShell* >::iterator _is = gwbasis.firstShell(); _is != gwbasis.lastShell(); _is++) {
                    AOShell* _shell = gwbasis.getShell(_is);
                    double decay = (*_shell->firstGaussian())->decay;
                    int _lmax    = _shell->getLmax();
                    int _size    = _shell->getNumFunc();
                    int _start  = _shell->getStartIndex();

                    const double pi = boost::math::constants::pi<double>();
                    double _factor = pow((2.0 *pi/decay),0.75);
                    vector<double> chi( _size, 0.0 );
                    chi[0] = _factor;

                    // some block from the fortran code that I'm not sure we need 
                    /*
                                  if ( lmax .ge. 0 ) then    
                      if(lmax .ge. 2 ) then
                       chi(10)= 6.d0*factor/sqrt(15.d0)   
                       if( lmax .ge. 4) then
                          fak = 0.25d0*factor*sqrt(beta_gwa) 
                          ! xxxx, yyyy, zzzz
                          chi(21) = fak*3.d0
                          chi(22) = chi(21)
                          chi(23) = chi(21)
                           ! xxzz, yyzz, xxyy
                           chi(30) = fak
                          chi(31) = fak
                          chi(32) = fak
                       end if
                    endif
                   end if
                     
                     */

                    // loop over all functions in shell
                    for ( int _i_gw = 0; _i_gw < _size ; _i_gw++ ){
                        double _test = _temp( _i_gw + _start, _n_level   );
                        if ( _test > 0.0  ){
                            sc_plus += chi[ _i_gw ]* _test;
                        } else if ( _test < 0.0 ){
                            sc_minus -= chi[ _i_gw ]* _test;
                        }
                    } // end loop over functions in shell

                } // end loop over all shells

                if ( _m_level <= _Mmn_RPA.get_mmax() && _n_level >= _Mmn_RPA.get_nmin()  ){
                    
                    double target = sqrt( sc_plus * sc_minus );
                    sc_plus  = target / sc_plus;
                    sc_minus = target / sc_minus;

                    // loop over gwbasis shells
                    for (vector< AOShell* >::iterator _is = gwbasis.firstShell(); _is != gwbasis.lastShell(); _is++) {
                        AOShell* _shell = gwbasis.getShell(_is);
                        int _size    = _shell->getNumFunc();
                        int _start  = _shell->getStartIndex();
                        vector<double> chi( _size, 0.0 );
                        chi[0] = 1.0;
                        // loop over all functions in shell
                        for ( int _i_gw = 0; _i_gw < _size ; _i_gw++ ){
                            double _test = _temp( _i_gw + _start, _n_level   );
                            if ( _test > 0.0 && std::abs( chi[_i_gw] ) > 1.e-10 ){
                               _temp( _i_gw + _start, _n_level   ) = _temp( _i_gw + _start, _n_level   ) * sc_plus;
                            } else if ( _test < 0.0 && std::abs( chi[_i_gw] ) > 1.e-10  ){
                               _temp( _i_gw + _start, _n_level   ) = _temp( _i_gw + _start, _n_level   ) * sc_minus;
                            }
                        } // end loop over functions in shell
                    } // end loop over all shells
                    
                }                
                
            }// loop n-levels

            // multiply _temp with overlap
            ub::matrix<double> _temp2 = ub::prod( gwoverlap._aomatrix , _temp );
            // copy to _Mmn_RPA
            _Mmn_RPA[ _m_level ] = ub::project( _temp2, ub::range(0, gwbasis._AOBasisSize) , ub::range(_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin()  , _Mmn_RPA.get_nmax() - _Mmn_full.get_nmin() +1 ));
            
            
        }// loop m-levels
        
    } // end RPA_prepare_threecenters


    }
    
 
};
