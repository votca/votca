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

#include "votca/ctp/elements.h"

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
            _numofelectrons =0;
            _mixingparameter=0.8;
            /*TEST for mkl and stuff
            ub::matrix<double> a=ub::zero_matrix<double>(3,3);
            ub::matrix<double> b=ub::zero_matrix<double>(3,3);
            ub::vector<double> E=ub::zero_vector<double>(3);
            ub::matrix<double> v=ub::zero_matrix<double>(3,3);
            a(0,0)=0;
            a(0,1)=2;
            a(0,2)=4;
            a(1,0)=2;
            a(2,0)=4;
            a(1,1)=4;
            a(2,1)=6;
            a(1,2)=6;
            a(2,2)=8;
            b(0,0)=3;
            b(0,1)=1;
            b(0,2)=1;
            b(1,0)=1;
            b(2,0)=1;
            b(1,1)=3;
            b(2,1)=1;
            b(1,2)=1;
            b(2,2)=3;
            linalg_eigenvalues_general(a,b,E,v);
            for (int i=0;i<3;i++){
                cout << "E("<< i << ")=" << E(i)<<endl;
                for (int j=0;j<3;j++){
                cout <<"v("<< j << ":"<< i <<")= " <<v(j,i)<< endl;
              //  cout <<"a("<< j << ":"<< i <<")= " <<a(j,i)<< endl;
                }}
            exit(0);
            */
            
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
            AOBasis* basis = &_dftbasis;
            int numofiterations=1;
           
	    /**** PREPARATION (atoms, basis sets, numerical integrations) ****/
	    Prepare( _orbitals );
            /**** END OF PREPARATION ****/
            
	    /**** Density-independent matrices ****/
	    SetupInvariantMatrices();
            
            int size4c=_dftbasis.AOBasisSize();
            fourdim fourcenter(boost::extents[size4c][size4c][size4c][size4c]);
            for ( int i = 0; i < size4c; i++ ){
                for ( int j = 0; j < size4c; j++ ){
                    for ( int k = 0; k < size4c; k++ ){
                        for ( int l = 0; l < size4c; l++ ){
                            
                            fourcenter[i][j][k][l] = 0.0;
                            
                        }
                    }
                }
            }
            
            LOG(logDEBUG, *_pLog) << TimeStamp() << "Init 4c "<< flush;
            
                         ifstream in;
            string _4cfile = "test/4cints";             
                       in.open(_4cfile.c_str(), ios::in);
                        if (!in) throw runtime_error(string("Error reading 4cs from: ")
                        + _4cfile);
                       string ao;
                       int iin;
                       int jin;
                       int kin;
                       int lin;
                       double cval;
                       int id = 0;
                       while (in.good()) { // keep reading until end-of-file
                           
                           id++;
                           //cout << "reading entry " << id << endl; 
                           
                           in >> ao;
                           in >> iin;
                           in >> jin;
                           in >> kin;
                           in >> lin;
                           in >> cval;
                              if (in.eof()) break;
                           // put to 4c
                           fourcenter[iin-1][jin-1][kin-1][lin-1] = cval;
                           
                           
                           
                       }
                       in.close();
                       
                       LOG(logDEBUG, *_pLog) << TimeStamp() << "Read 4cs from file "<< flush;
                       
                       
            /*    double x, y, z;
                int natoms, id;
                string label, type;
                vec pos;

                LOG(logDEBUG,_log) << " Reading molecular coordinates from " << _xyzfile << flush;
                in.open(_xyzfile.c_str(), ios::in);
                if (!in) throw runtime_error(string("Error reading coordinates from: ")
                        + _xyzfile);

                id = 1;
                while (in.good()) { // keep reading until end-of-file
                    in >> type;
                    in >> x;
                    in >> y;
                    in >> z;
                    if (in.eof()) break;
                    // cout << type << ":" << x << ":" << y << ":" << z << endl;

                    // creating atoms and adding them to the molecule
                    Atom *pAtom = new Atom(id++, type);
                    vec position(x / 10, y / 10, z / 10); // xyz has Angstrom, votca stores nm
                    pAtom->setPos(position);
                    pAtom->setQMPart(id, position);
                    pAtom->setElement(type);
                    _segment.AddAtom(pAtom);

                }
                in.close();

            */
            
            // _dftAOkinetic.Print("TMAT");
            //exit(0);
            
            /**** Initial guess = one-electron Hamiltonian without interactions ****/
            ub::vector<double>& MOEnergies=_orbitals->MOEnergies();
            ub::matrix<double>& MOCoeff=_orbitals->MOCoefficients();

            /**** Construct initial density  ****/
            
	               ub::matrix<double> H0=_dftAOkinetic._aomatrix+_dftAOESP._nuclearpotential;

		       //linalg_eigenvalues_general( H0,_dftAOoverlap._aomatrix, MOEnergies, MOCoeff);
            
		       /* double totinit=0;
         cout << endl;
                for (int i=0;i<(_numofelectrons/2);i++){
                    cout << MOEnergies(i) << " eigenwert " << i << endl;
                    totinit+=2*MOEnergies(i);
                }
                 LOG(logDEBUG, *_pLog) << TimeStamp() << " Total KS orbital Energy "<<totinit<<flush;
                 
		       */
            /*
            // cout << MOEnergies[0] << " " << MOEnergies[_numofelectrons/2] << endl;
            cout << "\n";
            for ( int i =0; i < _dftbasis.AOBasisSize(); i++){
            ub::matrix_range< ub::matrix<double> > _MO = ub::subrange(MOCoeff, 0, _dftbasis.AOBasisSize(), i, i+1 );
            
            ub::matrix<double> _temp = ub::prod(_dftAOoverlap._aomatrix, _MO);
            for ( int j =0; j < _dftbasis.AOBasisSize(); j++){
                ub::matrix_range< ub::matrix<double> > _MO2 = ub::subrange(MOCoeff, 0, _dftbasis.AOBasisSize(), j, j+1 );
            ub::matrix<double> norm = ub::prod(ub::trans(_MO2), _temp);
            if (std::abs(norm(0,0)) > 1e-5 ){
            cout << " Norm of MO " << i << " " << j << " " << norm(0,0) << endl;
            }}}
            exit(0);
            
            */
            
            
            
            
		       ub::matrix<double> initMOCoeff= ub::trans(_orbitals->MOCoefficients());
            
	    DensityMatrixGroundState( initMOCoeff, _numofelectrons/2 ) ;
	    cout << endl;
	 /*   for (int i =0 ; i < _dftAOdmat.size1(); i++){
	      for (int j =0 ; j < _dftAOdmat.size2(); j++){
		cout << "DMAT" << " " << i+1 << " " << j+1 << " " << _dftAOdmat(i,j) << endl;
	    }}
*/



            LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup Initial Guess "<< flush;
           LOG(logDEBUG, *_pLog) << TimeStamp() << " Num of electrons "<< _gridIntegration.IntegrateDensity(_dftAOdmat, basis) << flush;
	   
            for (int i=0;i<numofiterations;i++){
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Iteration "<< i+1 <<" of "<<numofiterations<< flush;


                // _ERIs.CalculateERIs(_dftAOdmat);
                // LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2()<< flush<<flush;


                /* ERI from 4cs*/
                ub::matrix<double> ERI4c = ub::zero_matrix<double>(size4c,size4c);
                for ( int ii=0; ii< size4c; ii++){
                    for ( int jj=0; jj< size4c; jj++){
                        for ( int kk=0; kk< size4c; kk++){
                            for ( int ll=0; ll< size4c; ll++){
                                
                                
                                ERI4c(ii,jj) += _dftAOdmat(kk,ll)*fourcenter[ii][jj][kk][ll];
                                
                                
                                
                            }
                        }
                    }
                }
                
                
		ub::matrix<double> VXC=_gridIntegration.IntegrateVXC_Atomblock(_dftAOdmat,  basis);
                //ub::matrix<double> H=H0+_ERIs.getERIs()+VXC;

                ub::matrix<double> H=H0+ERI4c+VXC;
                
                /*
		for ( int iout=0; iout<_dftAOdmat.size1();iout++){
		for ( int jout=0; jout<_dftAOdmat.size1();jout++){

		  cout.precision(10);
		  cout << "HAM " << std::setprecision(10)  << iout+1 << " " << jout+1 << " " << _dftAOkinetic._aomatrix(iout,jout) << " " << _dftAOESP._nuclearpotential(iout,jout) << " " << _ERIs.getERIs()(iout,jout) << " " << VXC(iout,jout)  << " " << _dftAOoverlap._aomatrix(iout,jout) << endl;

		}
		}
		//exit(0);
             */
		/*                for (int j=0;j<_ERIs.getSize2();j++){
                        for (int i=0;i<_ERIs.getSize1();i++){
                cout << "_ERIs ("<< i <<":"<< j<<")="<<_ERIs.getERIs()(i,j)<<endl;
		}}*/
            
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Vxc matrix "<<flush;
                linalg_eigenvalues_general( H,_dftAOoverlap._aomatrix, MOEnergies, MOCoeff);
                
                
               /* 		for ( int iout=0; iout<_dftAOdmat.size1();iout++){
		for ( int jout=0; jout<_dftAOdmat.size1();jout++){

		  cout.precision(10);
		  cout << "HAM2 " << std::setprecision(10)  << iout+1 << " " << jout+1 << " " << _dftAOkinetic._aomatrix(iout,jout) << " " << _dftAOESP._nuclearpotential(iout,jout) << " " << _ERIs.getERIs()(iout,jout) << " " << VXC(iout,jout)  << " " << _dftAOoverlap._aomatrix(iout,jout) << endl;

		}
		}
		exit(0);
*/                
                
                double totenergy=0;
                cout << endl;

                for (int i=0;i<_numofelectrons;i++){
                    cout << MOEnergies(i) << " eigenwert " << i << endl;
                    totenergy+=2*MOEnergies(i);
                }
                cout << " GAP " << MOEnergies(7)-MOEnergies(6) << endl;
                
                 LOG(logDEBUG, *_pLog) << TimeStamp() << " Total KS orbital Energy "<<totenergy<<flush;
                totenergy+=_gridIntegration.getTotEcontribution()-0.5*_ERIs.getERIsenergy();

                LOG(logDEBUG, *_pLog) << TimeStamp() << " Exc contribution "<<_gridIntegration.getTotEcontribution()<<flush;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " E_H contribution "<<-0.5*_ERIs.getERIsenergy()<<flush;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Total Energy "<<totenergy<<flush;
                
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Solved general eigenproblem "<<flush;




                EvolveDensityMatrix( MOCoeff, _numofelectrons/2 ) ;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Num of electrons "<< _gridIntegration.IntegrateDensity(_dftAOdmat, basis) << flush;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Updated Density Matrix "<<flush;
            }
          
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


            _dftAOkinetic.Initialize(_dftbasis.AOBasisSize());
            _dftAOkinetic.Fill(&_dftbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Kinetic energy matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
            
            _dftAOESP.Initialize(_dftbasis.AOBasisSize());
            _dftAOESP.Fillnucpotential(&_dftbasis, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT nuclear potential matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
            
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
            
            
            // prepare invariant part of electron repulsion integrals
            _ERIs.Initialize(_dftbasis, _auxbasis, _auxAOoverlap, _auxAOcoulomb);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup invariant parts of Electron Repulsion integrals " << flush;

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
	
           Elements _elements; 
            //set number of electrons and such
           _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
           
           for (int i=0;i<_atoms.size();i++){
               _numofelectrons+=_elements.getNucCrg(_atoms[i]->type);
           }
           LOG(logDEBUG, *_pLog) << TimeStamp() << " Total number of electrons: " << _numofelectrons << flush;
           _orbitals->setNumberOfElectrons(_numofelectrons);
           _orbitals->setNumberOfLevels(_numofelectrons/2,_dftbasis.AOBasisSize()-_numofelectrons/2);
           
           

            // _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
            
      }
      
      void DFTENGINE::DensityMatrixGroundState( ub::matrix<double>& _MOs, int occulevels ) {
     
     // first fill Density matrix, if required
    //  if ( _dmatGS.size1() != _basis_set_size ) {
          int size=max(_MOs.size1(),_MOs.size2());
          // cout << "Size " << size << " occ levels " << occulevels << endl;
        _dftAOdmat = ub::zero_matrix<double>(size,size);
        for ( int _i=0; _i < size; _i++ ){
            for ( int _j=0; _j < size; _j++ ){
                for ( int _level=0; _level < occulevels ; _level++ ){
                 
                    //_dftAOdmat(_i,_j) += 2.0 * _MOs( _level , _i ) * _MOs( _level , _j );
                    _dftAOdmat(_i,_j) += 2.0 * _MOs(  _i , _level) * _MOs(  _j, _level );
                 
                }
            }
         }
     //}
    }
      
      
      void DFTENGINE::EvolveDensityMatrix(ub::matrix<double>& MOCoeff, int occulevels){
          
      ub::matrix<double> dftdmat_old=_dftAOdmat;
      DensityMatrixGroundState(MOCoeff, occulevels);
      _dftAOdmat=_mixingparameter*_dftAOdmat+(1.0-_mixingparameter)*dftdmat_old;
      
      /*DIIS or mixing can be implemented here*/
      
      }
      

      
    
    }
};
