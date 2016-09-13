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
            
	    // exchange and correlation as in libXC
        
            _xc_functional_name = options->get(key + ".xc_functional").as<string> ();
        
	   
            _numofelectrons =0;
            
	    _max_iter = options->get(key + ".max_iterations").as<int>();
            
            
            if ( options->exists(key+".convergence")) {
                
            if ( options->exists(key+".convergence.energy")) {
                  _Econverged=options->get(key+".convergence.energy").as<double>();
               }
               else{
                    _Econverged=1e-7;
               }    
            
                
             if ( options->exists(key+".convergence.method")) {
               string method= options->get(key+".method").as<string>();
               if (method=="DIIS"){
                   _usediis=true;
               }
               else if(method=="mixing"){
                   _usediis=false;
               }
               else{
                   cout<< "WARNING method not known" <<endl;
                   _usediis=false;
               }
             }
             else{
                 _usediis=true;
             }
               if ( options->exists(key+".convergence.mixing")) {
                   _mixingparameter=options->get(key+".convergence.mixing").as<double>();
               }
               else{
                    _mixingparameter=0.5;
               }
            
            if ( options->exists(key+".convergence.maxout")) {
                    _maxout=options->get(key+".convergence.maxout").as<bool>();
               }
               else{
                     _maxout=false;
               }
             
              if ( options->exists(key+".convergence.diislength")) {
                   _histlength=options->get(key+".convergence.diislength").as<int>();
               }
               
               else{
                    _histlength=10;         
               }
               if ( options->exists(key+".convergence.usefockmat")) {
                   _usefmat=options->get(key+".convergence.usefockmat").as<bool>();
               }
               else{
                    _usefmat=false;
               } 
        }
        if(!_usediis){
                _histlength=1;
                _maxout=false;
                
                
            }
        else{
                _maxerrorindex=0;
                _maxerror=0.0;
            }
            
            return;
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
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Molecule Coordinates [A] "  << flush;
            for(unsigned i=0;i<_atoms.size();i++){
                LOG(logDEBUG, *_pLog) << "\t\t "<< _atoms[i]->type<<" "<<_atoms[i]->x<<" "<<_atoms[i]->y<<" "<<_atoms[i]->z<<" "<<flush;
            }
            
       
            
           
	    /**** PREPARATION (atoms, basis sets, numerical integrations) ****/
	    Prepare( _orbitals );
            /**** END OF PREPARATION ****/
            
	    /**** Density-independent matrices ****/
            SetupInvariantMatrices();
            
            
            
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


                _ERIs.CalculateERIs(_dftAOdmat, _AuxAOcoulomb_inv);
                LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2()<< flush<<flush;


		ub::matrix<double> VXC=_gridIntegration.IntegrateVXC_Atomblock(_dftAOdmat,  &_dftbasis,_xc_functional_name);
         
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
                    LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculation has converged up to "<<_Econverged<<". after "<< _this_iter<<" iterations." <<flush;
                    break;
                }
                else{
                    energyold=totenergy;
                }


                EvolveDensityMatrix(_orbitals,H ) ;
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

            {
	    // DFT AOOverlap matrix
	    _dftAOoverlap.Initialize(_dftbasis.AOBasisSize());
            _dftAOoverlap.Fill(&_dftbasis);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftAOoverlap.Dimension() << flush;

	    // check DFT basis for linear dependence
            linalg_eigenvalues(_dftAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            
            //Not brilliant but we need S-1/2 for DIIS and I do not want to calculate it each time
            ub::matrix<double> _diagS = ub::zero_matrix<double>(_eigenvectors.size1(),_eigenvectors.size1() );
            for ( unsigned _i =0; _i < _eigenvalues.size() ; _i++){

                _diagS(_i,_i) = 1.0/sqrt(_eigenvalues[_i]);
            }
            ub::matrix<double> _temp = ub::prod( _diagS, ub::trans(_eigenvectors));
             _Sminusonehalf = ub::prod(_eigenvectors,_temp );
            
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of DFT Overlap matrix : " << _eigenvalues[0] << flush;
      }
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
        
 
            { // this is just for info and not needed afterwards
            AOOverlap _auxAOoverlap;
            _auxAOoverlap.Initialize(_auxbasis.AOBasisSize());
            _auxAOoverlap.Fill(&_auxbasis);
            
            
            linalg_eigenvalues(_auxAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Smallest eigenvalue of AUX Overlap matrix : " << _eigenvalues[0] << flush;
             }
            
            
            
            
            
            // AUX AOcoulomb matrix
            //cout << " _auxAOcoulomb" <<endl;
            {
            AOCoulomb                           _auxAOcoulomb;
            _auxAOcoulomb.Initialize(_auxbasis.AOBasisSize());
            _auxAOcoulomb.Fill(&_auxbasis);
           
            //cout << _auxAOcoulomb._aomatrix<<endl;
            // _auxAOcoulomb.Print("COU");
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;
            _AuxAOcoulomb_inv=ub::zero_matrix<double>( _auxAOcoulomb.Dimension(), _auxAOcoulomb.Dimension()); 
            linalg_invert( _auxAOcoulomb.Matrix() , _AuxAOcoulomb_inv);
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Inverted AUX Coulomb matrix" << flush;
     
            }
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
      
      
      
      
      bool DFTENGINE::EvolveDensityMatrix(Orbitals* _orbitals,const ub::matrix<double>& H ){
          
      if(_errormatrixhist.size()>_histlength){
          delete _mathist[_maxerrorindex];
          delete _errormatrixhist[+_maxerrorindex];
              _mathist.erase(_mathist.begin()+_maxerrorindex);
              _errormatrixhist.erase(_errormatrixhist.begin()+_maxerrorindex);
          }
          
      
          
          
      bool has_converged=false;    
      _dftAOdmat=_orbitals->DensityMatrixGroundState(_orbitals->MOCoefficients());
      
      //Calculate errormatrix and orthogonalize
      ub::matrix<double>temp=ub::prod(H,_dftAOdmat);
      
      ub::matrix<double> errormatrix=ub::prod(temp,_dftAOoverlap._aomatrix);
      temp=ub::prod(_dftAOdmat,H);
      errormatrix-=ub::prod(_dftAOoverlap._aomatrix,temp);
     
      temp=ub::prod(errormatrix,_Sminusonehalf);
      errormatrix=ub::prod(ub::trans(_Sminusonehalf),temp);
      
      temp.resize(0,0);
      
      double max=linalg_getMax(errormatrix);
      ub::matrix<double>* old=new ub::matrix<double>;     
      if(_usefmat){
          *old=H;         
      }
      else{
          *old=_dftAOdmat;
      }
       _mathist.push_back(old);
      ub::matrix<double>* olderror=new ub::matrix<double>; 
      *olderror=errormatrix;
       _errormatrixhist.push_back(olderror);
       if(_maxout){
          double error=linalg_getMax(errormatrix);
          if (error>_maxerror){
              _maxerror=error;
              _maxerrorindex=_mathist.size();
          }
      } 
       
      if (max<0.1 && _this_iter>4 && _usediis){
          
          ub::matrix<double> B=ub::zero_matrix<double>(_mathist.size()+1);
          ub::vector<double> a=ub::zero_vector<double>(_mathist.size()+1);
          a(0)=-1;
          for (unsigned i=1;i<B.size1();i++){
              B(i,0)=-1;
          }
          #pragma omp parallel for
          for (unsigned i=1;i<B.size1();i++){
              for (unsigned j=1;j<=i;j++){
                  B(i,j)=linalg_traceofProd(*_errormatrixhist[i-1],ub::trans(*_errormatrixhist[j-1]));
                  if(i!=j){
                    B(j,i)=B(i,j);
                  }
              }
          }
          
          
          ub::vector<double> c;
          linalg_qrsolve(c, B, a);
          if(_usefmat){
                ub::matrix<double>H_guess=ub::zero_matrix<double>(_mathist[0]->size1()); 
                 for (unsigned i=0;i<_mathist.size();i++){  
                H_guess+=c(i+1)*(*_mathist[i]);
                 }
                _dftAOdmat=DmatfromFockmatrix(_orbitals,H_guess);
              }
          else{
             _dftAOdmat=ub::zero_matrix<double>(_mathist[0]->size1()); 
            for (unsigned i=0;i<_mathist.size();i++){  
            _dftAOdmat+=c(i+1)*(*_mathist[i]);
            }
          }
          
            
      }
      else if(_this_iter > 0 && _mathist.size()>0){
          if(_usefmat){
              
              ub::matrix<double>H_guess=_mixingparameter*H+(1.0-_mixingparameter)*(*(_mathist.back()));
              _dftAOdmat=DmatfromFockmatrix(_orbitals,H_guess);
          }
          else{
          _dftAOdmat=_mixingparameter*_dftAOdmat+(1.0-_mixingparameter)*(*(_mathist.back()));
          }
      }
      else if(max<_Econverged){
          has_converged=true;
      }
  
      return has_converged;
      }
      
      ub::matrix<double> DFTENGINE::DmatfromFockmatrix(Orbitals* _orbitals,const ub::matrix<double>&H){
          ub::vector<double> MOEnergies;
           ub::matrix<double> temp;
           linalg_eigenvalues_general( H,_dftAOoverlap._aomatrix, MOEnergies, temp);
           ub::matrix<double> MOCoeff=ub::trans(temp);
           
           
           return _orbitals->DensityMatrixGroundState(MOCoeff);
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
