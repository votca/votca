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


    


    }
    
 
};
