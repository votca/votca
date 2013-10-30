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

#include "gwbse.h"

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
        // GWBSE MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        void GWBSE::CleanUp() {

        }

        int GWBSE::NumFuncShell(string shell_type) {
            int _nbf;
            if (shell_type == "S") {
                _nbf = 1;
            } else if (shell_type == "P") {
                _nbf = 3;
            } else if (shell_type == "D") {
                _nbf = 5;
            } else if (shell_type == "SP") {
                _nbf = 4;
            } else if (shell_type == "SPD") {
                _nbf = 9;
            }
            return _nbf;
        }

        int GWBSE::OffsetFuncShell(string shell_type) {
            int _nbf;
            if (shell_type == "S") {
                _nbf = 0;
            } else if (shell_type == "P") {
                _nbf = 1;
            } else if (shell_type == "D") {
                _nbf = 4;
            } else if (shell_type == "SP") {
                _nbf = 0;
            } else if (shell_type == "SPD") {
                _nbf = 0;
            }
            return _nbf;
        }

        void GWBSE::Initialize(Property *options) {

            _maverick = (_nThreads == 1) ? true : false;

            /* obsolete string key = "options." + Identify();
            _jobfile = options->get(key + ".job_file").as<string>(); */

            string key = "options." + Identify() + ".job";
            _jobfile = options->get(key + ".file").as<string>();

            /*_do_input = false;
            _do_run = false;
            _do_parse = false;
            _do_trim = false;
    
            // conversion to GW 
            _do_convert = false;
            _do_gwbse_input = false;
            _do_gwbse_run = false;
            _do_gwbse_parse = false;
        

            string _package_xml = options->get(key+".package").as<string> ();

    
            string _tasks_string = options->get(key+".tasks").as<string> ();
            if (_tasks_string.find("input") != std::string::npos) _do_input = true;
            if (_tasks_string.find("run") != std::string::npos) _do_run = true;
            if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;    
            // GW-BSE tasks
            if (_tasks_string.find("convert") != std::string::npos) _do_convert = true;   
            if (_tasks_string.find("gwbse_setup") != std::string::npos) _do_gwbse_input = true;
            if (_tasks_string.find("gwbse_exec") != std::string::npos) _do_gwbse_run = true;    
            if (_tasks_string.find("gwbse_read") != std::string::npos) _do_gwbse_parse = true;
    
            string _store_string = options->get(key+".store").as<string> ();
            if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
            if (_store_string.find("qppert") != std::string::npos) _store_qppert = true;
            if (_store_string.find("qpdiag") != std::string::npos) _store_qpdiag = true;
            if (_store_string.find("singlets") != std::string::npos) _store_singlets = true;
            if (_store_string.find("triplets") != std::string::npos) _store_triplets = true;
    
            load_property_from_xml( _package_options, _package_xml.c_str() );    
            key = "package";
            _package = _package_options.get(key+".name").as<string> ();


   
            // only required, if GWBSE is to be run
            if ( _do_gwbse_input || _do_gwbse_run || _do_gwbse_parse ){
                key = "options." + Identify();
                string _gwpackage_xml = options->get(key+".gwpackage").as<string> ();
                load_property_from_xml( _gwpackage_options, _gwpackage_xml.c_str() );  
                key = "package";
                _gwpackage = _gwpackage_options.get(key+".name").as<string> ();
            }
    
    
            // register all QM packages (Gaussian, turbomole, nwchem))
            QMPackageFactory::RegisterAll(); */
            cout << "I'm supposed to initialize GWBSE";

        }

        Job::JobResult GWBSE::EvalJob(Topology *top, Job *job, QMThread *opThread) {

            cout << "Starting GW-BSE";
            Orbitals _orbitals;
            Job::JobResult jres = Job::JobResult();
            Property _job_input = job->getInput();
            list<Property*> lSegments = _job_input.Select("segment");

            vector < Segment* > segments;
            int segId = lSegments.front()->getAttribute<int>("id");
            string segType = lSegments.front()->getAttribute<string>("type");

            Segment *seg = top->getSegment(segId);
            assert(seg->Name() == segType);
            segments.push_back(seg);

            Logger* pLog = opThread->getLogger();
            LOG(logINFO, *pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush;

            // load the DFT data 
            string orb_file = (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
            string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());
            string edft_work_dir = "OR_FILES";
            string DIR = edft_work_dir + "/molecules_gwbse/" + frame_dir;
            std::ifstream ifs((DIR + "/" + orb_file).c_str());
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loading DFT data from " << DIR << "/" << orb_file << flush;
            boost::archive::binary_iarchive ia(ifs);
            ia >> _orbitals;
            ifs.close();
            string _dft_package = _orbitals.getQMpackage();
            LOG(logDEBUG, *pLog) << TimeStamp() << " DFT data was created by " << _dft_package << flush;

            // reorder DFT data, load DFT basis set
            BasisSet dftbs;
            string dftbasis_name("ubecppol");

            AOBasis dftbasis;

            dftbs.LoadBasisSet(dftbasis_name);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loaded DFT Basis Set " << dftbasis_name << flush;

            dftbasis.AOBasisFill(&dftbs, segments);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled DFT Basis of size " << dftbasis._AOBasisSize << flush;

            // do the reordering depending on the QM package used to obtain the DFT data
            ub::matrix<double> _dft_orbitals = *_orbitals.getOrbitals();
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
            LOG(logDEBUG, *pLog) << TimeStamp() << " Converted DFT orbital coefficient order " << flush;

            // setting up ao_overlap_matrix
            list<string> elements;
            BasisSet gwbs;
            string gwbasis_name("gwdefault");

            AOBasis gwbasis;
            bool PPM_symmetric = true; // only PPM supported


            gwbs.LoadBasisSet(gwbasis_name);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loaded GW Basis Set " << gwbasis_name << flush;

            gwbasis.AOBasisFill(&gwbs, segments);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Basis of size " << gwbasis._AOBasisSize << flush;

            // get overlap matrix as AOOverlap
            AOOverlap _gwoverlap;
            // initialize overlap matrix
            _gwoverlap.Initialize(gwbasis._AOBasisSize);
            // Fill overlap
            _gwoverlap.Fill(&gwbasis);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Overlap matrix of dimension: " << _gwoverlap._aomatrix.size1() << flush;
            // _aooverlap.Print( "S" );



            // printing some debug info
            // _gwcoulomb.PrintIndexToFunction( &aobasis );

            // check eigenvalues of overlap matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            _eigenvalues.resize(_gwoverlap._aomatrix.size1());
            _eigenvectors.resize(_gwoverlap._aomatrix.size1(), _gwoverlap._aomatrix.size1());
            linalg_eigenvalues(_gwoverlap._aomatrix, _eigenvalues, _eigenvectors);
            // cout << _eigenvalues << endl;
            sort(_eigenvalues.begin(), _eigenvalues.end());
            LOG(logDEBUG, *pLog) << TimeStamp() << " Smallest eigenvalue of GW Overlap matrix : " << _eigenvalues[0] << flush;



            // get Coulomb matrix as AOCoulomb
            AOCoulomb _gwcoulomb;
            // initialize Coulomb matrix
            _gwcoulomb.Initialize(gwbasis._AOBasisSize);
            // Fill Coulomb matrix
            _gwcoulomb.Fill(&gwbasis);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Coulomb matrix of dimension: " << _gwcoulomb._aomatrix.size1() << flush;
            //_gwcoulomb.Print( "COU_in " );


            // PPM is symmetric, so we need to get the sqrt of the Coulomb matrix
            AOOverlap _gwoverlap_inverse;               // will also be needed in PPM itself
            AOOverlap _gwoverlap_cholesky_inverse;      // will also be needed in PPM itself
            if (PPM_symmetric) {
                _gwoverlap_inverse.Initialize( gwbasis._AOBasisSize);
                _gwcoulomb.Symmetrize(  _gwoverlap , gwbasis, _gwoverlap_inverse , _gwoverlap_cholesky_inverse );
            }
            LOG(logDEBUG, *pLog) << TimeStamp() << " Prepared GW Coulomb matrix for symmetric PPM " << flush;

            // calculate 3-center integrals,  convoluted with DFT eigenvectors

            // --- prepare a vector (gwdacay) of matrices (orbitals, orbitals) as container => M_mn
            mmin = 1; // lowest index occ 
            mmax = 2 * _orbitals.getNumberOfElectrons();
            nmin = 1;
            nmax = _orbitals.getNumberOfLevels();
            maxf = gwbasis.getMaxFunctions(); // maximum number of functions per shell in basis set
            mtotal = mmax - mmin + 1;
            ntotal = nmax - nmin + 1;


            // prepare 3-center integral object
            TCMatrix _Mmn;
            _Mmn.Initialize(gwbasis._AOBasisSize, mmin, mmax, nmin, nmax);
            _Mmn.Fill(gwbasis, dftbasis, _dft_orbitals);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Calculated Mmn_beta (3-center-overlap x orbitals)  " << flush;
            //_Mmn.Print( "Mmn " );
            
            // for use in RPA, make a copy of _Mmn with dimensions (1:HOMO)(gwabasissize,LUMO:nmax)
            TCMatrix _Mmn_RPA;
            _Mmn_RPA.Initialize(gwbasis._AOBasisSize, mmin, _orbitals.getNumberOfElectrons() , _orbitals.getNumberOfElectrons() +1 , nmax);
            RPA_prepare_threecenters( _Mmn_RPA, _Mmn, gwbasis, _gwoverlap, _gwoverlap_inverse );
            

            // TODO: now, we can get rid of _gwoverlap_inverse
            // make _Mmn_RPA symmetric for use in RPA
            _Mmn_RPA.Symmetrize( _gwcoulomb._aomatrix  );
            LOG(logDEBUG, *pLog) << TimeStamp() << " Prepared Mmn_beta for RPA  " << flush;
            // _Mmn_RPA.Print( "Mmn_RPA" );
            
            // make _Mmn symmetric for use in self-energy calculation
            _Mmn.Symmetrize( _gwcoulomb._aomatrix  );
            LOG(logDEBUG, *pLog) << TimeStamp() << " Prepared Mmn_beta for self-energy  " << flush;


            // some parameters that need to be prepared by options parsing, here fixed for testing
            _shift = 0.3; // in Rydberg
            _screening_freq = ub::zero_matrix<double>(2,2); // two frequencies
            //first one
            _screening_freq(0,0) = 0.0; // real part
            _screening_freq(0,1) = 0.0; // imaginary part
            //second one
            _screening_freq(1,0) = 0.0; // real part
            _screening_freq(1,1) = 1.0; // imaginary part

            ub::vector<double> _dft_energies = 2.0*(*_orbitals.getEnergies()); // getEnergies -> Hartree, we want Ryd
            
            // one entry to epsilon for each frequency
            _epsilon.resize( _screening_freq.size1() );
            
            // for symmetric PPM, we can initialize _epsilon with the overlap matrix!
            for ( int _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
                _epsilon( _i_freq ) = _gwoverlap._aomatrix ;
            }
            // TODO: now, we can get rid of _gwoverlap
            
            // determine epsilon from RPA
            RPA_calculate_epsilon( _Mmn_RPA, _screening_freq, _shift, _dft_energies );
            LOG(logDEBUG, *pLog) << TimeStamp() << " Calculated epsilon via RPA  " << flush;
  
            // construct PPM parameters
            PPM_construct_parameters( _gwoverlap_cholesky_inverse._aomatrix );
            
            
            
            LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
 
            Property _job_summary;
            Property *_output_summary = &_job_summary.add("output","");
            Property *_segment_summary = &_output_summary->add("segment","");
            string segName = seg->getName();
            segId = seg->getId();
            _segment_summary->setAttribute("id", segId);
            _segment_summary->setAttribute("type", segName);
            // output of the JOB 
            jres.setOutput( _job_summary );
            jres.setStatus(Job::COMPLETE);

            // dump the LOG
            cout << *pLog;
            return jres;
        }

        
        
        void GWBSE::PPM_construct_parameters(  ub::matrix<double>& _overlap_cholesky_inverse ){
            
            // multiply with L-1^t from the right
            ub::matrix<double> _overlap_cholesky_inverse_transposed = ub::trans( _overlap_cholesky_inverse );
            ub::matrix<double> _temp = ub::prod( _epsilon(0) , _overlap_cholesky_inverse_transposed );
            // multiply with L-1 from the left
            _temp = ub::prod( _overlap_cholesky_inverse, _temp );
            
            // get eigenvalues and eigenvectors of this matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            _eigenvalues.resize(_temp.size1());
            _eigenvectors.resize(_temp.size1(), _temp.size1());
            linalg_eigenvalues(_temp, _eigenvalues, _eigenvectors);
            // cout << _eigenvalues << endl;
            //sort(_eigenvalues.begin(), _eigenvalues.end());
            // cout << "min EV: " << _eigenvalues(0) << " max EV: " << _eigenvalues(_temp.size1() -1 ) << endl;
            
            
            // multiply eigenvectors with overlap_cholesky_inverse_transpose and store as eigenvalues of epsilon
            _ppm_phi = ub::prod( _overlap_cholesky_inverse_transposed , _eigenvectors ); 
            



            // store PPM weights from eigenvalues
            _ppm_weight.resize( _eigenvalues.size() );
            for ( int _i = 0 ; _i <  _eigenvalues.size(); _i++   ){
                _ppm_weight(_i) = 1.0 - 1.0/_eigenvalues(_i);
		// cout << " PPM weight " << _eigenvalues(_i) << " : " << _ppm_weight(_i) << endl; 
            }
            



            // determine PPM frequencies
            _ppm_freq.resize( _eigenvalues.size() );
            // a) phi^t * epsilon(1) * phi 
            _temp = ub::prod( ub::trans( _ppm_phi ) , _epsilon(1) );
            _eigenvectors  = ub::prod( _temp ,  _ppm_phi  );


	    // before inversion


	    /*	    for ( int i =0 ; i < _eigenvalues.size() ; i++){
	      for ( int j =0 ; j < _eigenvalues.size() ; j++){
	  
		cout << " PPM_before " << i << ":" << j << "  " << _eigenvectors(i,j) << endl;

	      }

	      } */ 




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
            
                   
            
        }
        
        
        
        
        
        

        void GWBSE::RPA_calculate_epsilon(TCMatrix& _Mmn_RPA, ub::matrix<double> _screening_freq, double _shift, ub::vector<double>& _dft_energies){
            
            int _size = _Mmn_RPA.matrix()(0).size1(); // size of gwbasis
            
            // loop over frequencies
            for ( int _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
                // initialize epsilon for this frequency
                // _epsilon ( _i_freq ) = ub::zero_matrix<double>(_size, _size);
                
                // loop over occupied bands -> vector index of _Mmn_RPA
                for ( int _m_band = 0; _m_band < _Mmn_RPA.get_mtot() ; _m_band++ ){
                    int index_m = _Mmn_RPA.get_mmin();

                    // a temporary matrix, that will get filled in empty bands loop
                    ub::matrix<double> _temp = ub::zero_matrix<double>( _Mmn_RPA.get_ntot(), _size );
                    
                        
                    // loop over empty bands
                    for ( int _n_band = 0 ; _n_band < _Mmn_RPA.get_ntot() ; _n_band++ ){
                        int index_n = _Mmn_RPA.get_nmin();
                        
                        
                        double _deltaE = _shift + _dft_energies( _n_band + index_n ) - _dft_energies( _m_band + index_m ); // get indices and units right!!!
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

                        for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                            _temp( _n_band , _i_gw ) = _energy_factor * _Mmn_RPA.matrix()( _m_band )( _i_gw , _n_band );
                        } // matrix size
                        
                    } // empty bands

                   // now multiply and add to epsilon
                   _epsilon( _i_freq ) += ub::prod( _Mmn_RPA.matrix()( _m_band ) , _temp  );

                } // occupied bands
                
            } // loop over frequencies
            
            
        }
        
        
   
    
    void GWBSE::RPA_prepare_threecenters( TCMatrix& _Mmn_RPA, TCMatrix& _Mmn_full, AOBasis& gwbasis, AOMatrix& gwoverlap, AOMatrix& gwoverlap_inverse     ){
        
        // cout << "blabla" << endl;
        
        
        // loop over m-bands in _Mmn_full
        // for ( int _m_band = 0; _m_band < _Mmn_full.matrix().size() ; _m_band++ ){
        // actually, only needed for size of _Mmn_RPA (until VBM)
        for ( int _m_band = 0; _m_band < _Mmn_RPA.matrix().size() ; _m_band++ ){
        
            ub::matrix<double> _temp = ub::prod( gwoverlap_inverse._aomatrix , _Mmn_full.matrix()( _m_band ) );

            // loop over n-bands in _Mmn_full 
            for ( int _n_band = 0; _n_band < _Mmn_full.get_ntot() ; _n_band++ ){

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
                        double _test = _temp( _i_gw + _start, _n_band   );
                        if ( _test > 0.0  ){
                            sc_plus += chi[ _i_gw ]* _test;
                        } else if ( _test < 0.0 ){
                            sc_minus -= chi[ _i_gw ]* _test;
                        }
                    } // end loop over functions in shell

                } // end loop over all shells

                if ( _m_band <= _Mmn_RPA.get_mmax() && _n_band >= _Mmn_RPA.get_nmin()  ){
                    
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
                            double _test = _temp( _i_gw + _start, _n_band   );
                            if ( _test > 0.0 && std::abs( chi[_i_gw] ) > 1.e-10 ){
                               _temp( _i_gw + _start, _n_band   ) = _temp( _i_gw + _start, _n_band   ) * sc_plus;
                            } else if ( _test < 0.0 && std::abs( chi[_i_gw] ) > 1.e-10  ){
                               _temp( _i_gw + _start, _n_band   ) = _temp( _i_gw + _start, _n_band   ) * sc_minus;
                            }
                        } // end loop over functions in shell
                    } // end loop over all shells
                    
                }                
                
            }// loop n-bands

            // multiply _temp with overlap
            ub::matrix<double> _temp2 = ub::prod( gwoverlap._aomatrix , _temp );
            // copy to _Mmn_RPA
            // range(start, stop) contains all indices i with start <= i < stop
            _Mmn_RPA.matrix()( _m_band ) = ub::project( _temp2, ub::range(0, gwbasis._AOBasisSize) , ub::range(_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin()  , _Mmn_RPA.get_nmax() - _Mmn_full.get_nmin() +1 ));
            
            
        }// loop m-bands
        
    } // end RPA_prepare_threecenters


    }
    
 
};
