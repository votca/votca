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

#include <votca/ctp/dftengine.h>

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

// #include <omp.h>



using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // DFTENGINE MEMBER FUNCTIONS        //
        // +++++++++++++++++++++++++++++ //

        void DFTENGINE::CleanUp() {

        }
        
        
        void DFTENGINE::Initialize(Property* options){

           // setting some defaults
           // _do_qp_diag = false;
            _openmp_threads = 0; // take all available
            
            string key =  Identify();
    
            // get OpenMP thread number
            _openmp_threads = options->get(key + ".openmp").as<int> ();
            
            // basis sets
	    _dftbasis_name = options->get(key + ".dftbasis").as<string>();
	    _auxbasis_name = options->get(key + ".auxbasis").as<string>();

	    // numerical integrations
	    _grid_name = options->get(key + ".integration_grid").as<string>();

	    // exchange and correlation as in libXC
	    _x_functional_name = options->get(key + ".exchange_functional").as<string>();
	    _c_functional_name = options->get(key + ".correlation_functional").as<string>();
            
        }
        
        
        
        
        
        /* 
         *    Density Functional theory implementation
         * 
         */
        
        bool DFTENGINE::Evaluate( Orbitals* _orbitals) {

            
            // set the parallelization 
#ifdef OMP
            if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads);
#endif
            _atoms = _orbitals->QMAtoms();


	    /**** PREPARATION (atoms, basis sets, numerical integrations) ****/
	    Prepare( _orbitals );
            /**** END OF PREPARATION ****/
            
	    /**** Density-independent matrices ****/
	    SetupInvariantMatrices();

            /**** Initial guess = one-electron Hamiltonian without interactions ****/
            
            /**** Construct initial density  ****/


       
            /* STUPID TEST
            
            
            // test numerical integration for AOoverlap
            ub::matrix<double> OLTEST = _numint.numAOoverlap(&ddftbasis);
            
            // test total number of electrons
            ub::matrix<double> _ddft_orbitals = *(_orbitals->getOrbitals()); 
            ddftbasis.ReorderMOs(_ddft_orbitals, "gaussian", "votca" );
            ub::matrix<double> &DMAT = _orbitals->DensityMatrixGroundState( _ddft_orbitals );
            double Nelectrons = _numint.IntegrateDensity(DMAT,&ddftbasis);
            cout << " Number of electrons: " << Nelectrons << endl;

            // test AOxcmatrix
            ub::matrix<double> AOXC = _numint.IntegrateVXC(DMAT,&ddftbasis); 
            for ( int i = 0 ; i < AOXC.size1(); i++ ){
                //for ( int j = 0 ; j < AOXC.size2(); j++ ){
                    
                    cout << i << " : " << i << " aoxc " << AOXC(i,i) << endl;
                    
               // }
                
            }


            
            // get integration grid as matrix (needed?)
            ub::matrix<double> _intgrid;
            _numint.getGridpoints( _intgrid );
            
            //_gwoverlap.Print("AOOL");
            

	    */
            
            
            return true;
        }



      // SETUP INVARIANT AOMatrices
      void DFTENGINE::SetupInvariantMatrices(){


	// local variables for checks
            // check eigenvalues of overlap matrix, if too small basis might have linear dependencies
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;


	    // DFT AOOverlap matrix
	    _dftAOoverlap.Initialize(_dftbasis.AOBasisSize());
            _dftAOoverlap.Fill(&_dftbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftAOoverlap.Dimension() << flush;

	    // check DFT basis for linear dependence
            linalg_eigenvalues(_dftAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of DFT Overlap matrix : " << _eigenvalues[0] << flush;

	    // AUX AOoverlap
            _auxAOoverlap.Initialize(_auxbasis.AOBasisSize());
            _auxAOoverlap.Fill(&_auxbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled AUX Overlap matrix of dimension: " << _auxAOoverlap.Dimension() << flush;
 
	    // check AUX basis for linear dependence
            linalg_eigenvalues(_auxAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of AUX Overlap matrix : " << _eigenvalues[0] << flush;

            
            // AUX AOcoulomb matrix
            _auxAOcoulomb.Initialize(_auxbasis.AOBasisSize());
            _auxAOcoulomb.Fill(&_auxbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;


      }




      // PREPARATION 
      void DFTENGINE::Prepare( Orbitals* _orbitals ){

            // load and fill DFT basis set
            _dftbasisset.LoadBasisSet(_dftbasis_name);
            _orbitals->setDFTbasis( _dftbasis_name );    
	    _dftbasis.AOBasisFill( &_dftbasisset, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

	    // load and fill AUX basis set
            _auxbasisset.LoadBasisSet(_auxbasis_name);
            //_orbitals->setDFTbasis( _dftbasis_name );
	    _auxbasis.AOBasisFill( &_auxbasisset, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << flush;
	    
	    // setup numerical integration grid
            _gridIntegration.GridSetup(_grid_name,&_dftbasisset,_atoms);
	    LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup numerical integration grid " << _grid_name << flush;
	

      }


    }
    
 
};
