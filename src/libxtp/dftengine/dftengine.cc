/*
 *            Copyright 2009-2016 The VOTCA Development Team
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
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/dftengine.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/threecenters.h>
// #include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

#include "votca/xtp/elements.h"

// #include <omp.h>



using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
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
             if ( options->exists(key+".ecp")) {
               _ecp_name = options->get(key+".ecp").as<string>();
               _with_ecp = true;
             } else {
                 _with_ecp = false;
             }
            
             if ( options->exists(key+".guess")) {
               _ecp_name = options->get(key+".guess").as<bool>();
               _with_guess = true;
             } else {
                 _with_guess = false;
             }
            
            
	    // numerical integrations
	    _grid_name = options->get(key + ".integration_grid").as<string>();
            _Econverged = options->get(key + ".convergence").as<double>();
	    // exchange and correlation as in libXC
        
            _xc_functional_name = options->get(key + ".xc_functional").as<string> ();
        
	   
            _numofelectrons =0;
            _mixingparameter = options->get(key + ".density_mixing").as<double>();
	    _max_iter = options->get(key + ".max_iterations").as<int>();
            
            
        }
        
        
        
        
        
        /* 
         *    Density Functional theory implementation
         * 
         */
        
        bool DFTENGINE::Evaluate( Orbitals* _orbitals) {

            
            // set the parallelization 
            #ifdef _OPENMP
            if ( _openmp_threads > 0 ) {
                omp_set_num_threads(_openmp_threads);
                LOG(logDEBUG, *_pLog) << TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
            }
            #endif

            _atoms = _orbitals->QMAtoms();
            LOG(logDEBUG, *_pLog) << TimeStamp() << "Molecule Coordinates [A] "  << flush;
            for(unsigned i=0;i<_atoms.size();i++){
                LOG(logDEBUG, *_pLog) << "\t\t"<< _atoms[i]->type<<" "<<_atoms[i]->x<<" "<<_atoms[i]->y<<" "<<_atoms[i]->z<<" "<<flush;
            }
            
            AOBasis* basis = &_dftbasis;
            
           
	    /**** PREPARATION (atoms, basis sets, numerical integrations) ****/
	    Prepare( _orbitals );
            /**** END OF PREPARATION ****/
            
	    /**** Density-independent matrices ****/
            SetupInvariantMatrices();
            ub::matrix<double> _inverse_Coulomb=ub::zero_matrix<double>( _auxAOcoulomb.Dimension(), _auxAOcoulomb.Dimension()); 
            linalg_invert( _auxAOcoulomb.Matrix() , _inverse_Coulomb);
            
            
            /**** Initial guess = one-electron Hamiltonian without interactions ****/
            ub::vector<double>& MOEnergies=_orbitals->MOEnergies();
            ub::matrix<double>& MOCoeff=_orbitals->MOCoefficients();

            /**** Construct initial density  ****/

            ub::matrix<double> H0 = _dftAOkinetic._aomatrix + _dftAOESP._nuclearpotential; 
            
            if(_with_ecp){
            H0+=_dftAOECP.Matrix();
            cout<< "WARNING ecps are not correctly sorted" <<endl;
            }
            
            // if we have a guess we do not need this. 
            if(!_with_guess){
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup Initial Guess "<< flush;
             // this temp is necessary because eigenvalues_general returns MO^T and not MO
            ub::matrix<double> temp;
            linalg_eigenvalues_general(H0, _dftAOoverlap._aomatrix, MOEnergies,temp);
            MOCoeff=ub::trans(temp);
            }
            double totinit = 0;

            for (int i = 0; i < (_numofelectrons / 2); i++) {
                //cout << MOEnergies(i) << " eigenwert " << i << endl;
                totinit += 2 * MOEnergies(i);
            }
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Total KS orbital Energy " << totinit << flush;
            
            NuclearRepulsion();
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Nuclear Repulsion Energy is " << E_nucnuc << flush;

         
            _dftAOdmat=_orbitals->DensityMatrixGroundState(MOCoeff);
	    

        
	    


           
           //LOG(logDEBUG, *_pLog) << TimeStamp() << " Num of electrons "<< _gridIntegration.IntegrateDensity_Atomblock(_dftAOdmat, basis) << flush;
	   
           double energyold=totinit;
           
            for ( _this_iter=0; _this_iter<_max_iter; _this_iter++){
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Iteration "<< _this_iter+1 <<" of "<< _max_iter << flush;


                _ERIs.CalculateERIs(_dftAOdmat, _inverse_Coulomb);
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2()<< flush<<flush;


		ub::matrix<double> VXC=_gridIntegration.IntegrateVXC_Atomblock(_dftAOdmat,  basis,_xc_functional_name);
         
       //cout << "ERIS"<<endl;
       //cout<<_ERIs.getERIs()<<endl;
                
                
                ub::matrix<double> H=H0+_ERIs.getERIs()+VXC;
             
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Vxc matrix "<<flush;
                
                // this temp is necessary because eigenvalues_general returns MO^T and not MO
                ub::matrix<double> temp;
                linalg_eigenvalues_general( H,_dftAOoverlap._aomatrix, MOEnergies, temp);
                MOCoeff=ub::trans(temp);
          
                
                double totenergy=E_nucnuc;
       

                for (int i=0;i<_numofelectrons;i++){
                    if ( i <= _numofelectrons/2-1) {
                        LOG(logDEBUG, *_pLog) <<"\t\t" << i <<  " occ " << MOEnergies(i)  << flush;
                        totenergy+=2*MOEnergies(i);
                    } else {
                        LOG(logDEBUG, *_pLog) <<"\t\t"<<   i <<  " vir " << MOEnergies(i)  << flush;
                        
                    }
                }
                LOG(logDEBUG, *_pLog) << "\t\tGAP " << MOEnergies(_numofelectrons/2)-MOEnergies(_numofelectrons/2-1) << flush;
                
                 LOG(logDEBUG, *_pLog) << TimeStamp() << " Total KS orbital Energy "<<totenergy<<flush;
                totenergy+=_gridIntegration.getTotEcontribution()-0.5*_ERIs.getERIsenergy();

                LOG(logDEBUG, *_pLog) << TimeStamp() << " Exc contribution "<<_gridIntegration.getTotEcontribution()<<flush;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " E_H contribution "<<0.5*_ERIs.getERIsenergy()<<flush;
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Total Energy "<<totenergy<<flush;
                
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Solved general eigenproblem "<<flush;
                if (std::abs(totenergy-energyold)< _Econverged){
                    LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculation has converged up to "<<_Econverged<<"."<<flush;
                    break;
                }
                else{
                    energyold=totenergy;
                }


                EvolveDensityMatrix(_orbitals ) ;
                //LOG(logDEBUG, *_pLog) << TimeStamp() << " Num of electrons "<< _gridIntegration.IntegrateDensity_Atomblock(_dftAOdmat, basis) << flush;
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
            _dftAOESP.Fillnucpotential(&_dftbasis, _atoms,_with_ecp);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT nuclear potential matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
            //_dftAOESP.Print("NUC");

            if (_with_ecp) {
                _dftAOECP.Initialize(_dftbasis.AOBasisSize());
                _dftAOECP.Fill(&_dftbasis, ub::zero_vector<double>(3), &_ecp);
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT ECP matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
                //_dftAOECP.Print("ECP");
                
                _dftAOESP._nuclearpotential += _dftAOECP.Matrix();
                
            }
            // exit(0);
            
	    // AUX AOoverlap
        
 
   
            _auxAOoverlap.Initialize(_auxbasis.AOBasisSize());
            _auxAOoverlap.Fill(&_auxbasis);
            
            
            linalg_eigenvalues(_auxAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of AUX Overlap matrix : " << _eigenvalues[0] << flush;
 
            
            
            
            
            
            // AUX AOcoulomb matrix
            //cout << " _auxAOcoulomb" <<endl;
            _auxAOcoulomb.Initialize(_auxbasis.AOBasisSize());
            _auxAOcoulomb.Fill(&_auxbasis);
            
            //cout << _auxAOcoulomb._aomatrix<<endl;
            // _auxAOcoulomb.Print("COU");
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;
            //exit(0);
       
            // prepare invariant part of electron repulsion integrals
      
            _ERIs.Initialize(_dftbasis, _auxbasis);
          
            
            
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup invariant parts of Electron Repulsion integrals " << flush;

      }




      // PREPARATION 
      void DFTENGINE::Prepare( Orbitals* _orbitals ){

            // load and fill DFT basis set
            _dftbasisset.LoadBasisSet(_dftbasis_name);
            
            if(_with_guess && _orbitals->getDFTbasis()!=_dftbasis_name){
                throw runtime_error("Basisset Name in guess orb file and in dftengine option file differ.");
            }
            _orbitals->setDFTbasis( _dftbasis_name );    
	    _dftbasis.AOBasisFill( &_dftbasisset, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

	    // load and fill AUX basis set
            _auxbasisset.LoadBasisSet(_auxbasis_name);
            //_orbitals->setDFTbasis( _dftbasis_name );
	    _auxbasis.AOBasisFill( &_auxbasisset, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << flush;

            if(_with_ecp){
            // load ECP (element-wise information) from xml file
            _ecpbasisset.LoadPseudopotentialSet("ecp");
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Loaded ECP library " << "ecp" << flush;

            // fill auxiliary ECP basis by going through all atoms
            _ecp.ECPFill(&_ecpbasisset, _atoms);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled ECP Basis of size " << _ecp._aoshells.size() << flush;
            }
            
	    // setup numerical integration grid
            _gridIntegration.GridSetup(_grid_name,&_dftbasisset,_atoms);
	    LOG(logDEBUG, *_pLog) << TimeStamp() << " Setup numerical integration grid " << _grid_name << " for vxc functional " << _xc_functional_name<< flush;
        
	
           Elements _elements; 
            //set number of electrons and such
           _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
           
           for (unsigned i=0;i<_atoms.size();i++){
               _numofelectrons+=_elements.getNucCrg(_atoms[i]->type);
           }

            // if ECP
            if (_with_ecp) {
                for (unsigned i = 0; i < _atoms.size(); i++) {
                    _numofelectrons-= _ecpbasisset.getElement(_atoms[i]->type)->getNcore() ;
                }
            }
           
           
           LOG(logDEBUG, *_pLog) << TimeStamp() << " Total number of electrons: " << _numofelectrons << flush;
           
           if(_with_guess){
               if(_orbitals->getNumberOfElectrons()!=_numofelectrons){
               throw runtime_error("Number of electron in guess orb file and in dftengine differ.");
               }
               if(_orbitals->getNumberOfLevels()!=_dftbasis.AOBasisSize()){
               throw runtime_error("Number of levels in guess orb file and in dftengine differ.");
               }
            }
           
           _orbitals->setNumberOfElectrons(_numofelectrons);
           _orbitals->setNumberOfLevels(_numofelectrons/2,_dftbasis.AOBasisSize()-_numofelectrons/2);
           
           

            // _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
            
      }
      
      
      
      
      void DFTENGINE::EvolveDensityMatrix(Orbitals* _orbitals ){
          
      ub::matrix<double> dftdmat_old=_dftAOdmat;
      _dftAOdmat=_orbitals->DensityMatrixGroundState(_orbitals->MOCoefficients());;
      if (_this_iter > 0) _dftAOdmat=_mixingparameter*_dftAOdmat+(1.0-_mixingparameter)*dftdmat_old;
      
      /*DIIS or mixing can be implemented here*/
      
      }
      
      void DFTENGINE::NuclearRepulsion(){
          Elements element;
          E_nucnuc=0.0;
          
          std::vector<double> charge;
          for(unsigned i=0;i<_atoms.size();i++){
              if(_with_ecp){
                  charge.push_back(element.getNucCrgECP(_atoms[i]->type));
              }
              else{
                  charge.push_back(element.getNucCrg(_atoms[i]->type));
              }
          }      
              
          for(unsigned i=0;i<_atoms.size();i++){
              vec r1=vec(_atoms[i]->x*tools::conv::ang2bohr,_atoms[i]->y*tools::conv::ang2bohr,_atoms[i]->z*tools::conv::ang2bohr);
              double charge1=charge[i];
              for(unsigned j=0;j<i;j++){
                  vec r2=vec(_atoms[j]->x*tools::conv::ang2bohr,_atoms[j]->y*tools::conv::ang2bohr,_atoms[j]->z*tools::conv::ang2bohr);
                  double charge2=charge[j];
                  E_nucnuc+=charge1*charge2/(abs(r1-r2));
              }
          }

          return;
      }
      
      

      
    
    }
};
