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
#include <votca/tools/linalg.h>
#include "votca/xtp/aobasis.h"
#include <votca/xtp/dftengine.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/threecenters.h>

#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

#include <votca/xtp/elements.h>
#include <votca/xtp/diis.h>

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
            
            if ( options->exists(key+".auxbasis")) {
                _auxbasis_name = options->get(key + ".auxbasis").as<string>();
                _with_RI=true;
               
             } else {
                 
                 _with_RI=false;
             }
            
            if ( options->exists(key+".fourcmethod")) {
                _4cmethod = options->get(key + ".fourcmethod").as<string>();
                //direct or ram
               
             } else {
                 
                _4cmethod="ram";
             }
	   
             if ( options->exists(key+".ecp")) {
               _ecp_name = options->get(key+".ecp").as<string>();
               _with_ecp = true;
             } else {
                 _with_ecp = false;
             }
            
             if ( options->exists(key+".read_guess")) {
               _with_guess = options->get(key+".read_guess").as<bool>();
             } else {
                 _with_guess = false;
             }
             if ( options->exists(key+".initial_guess")) {
               _initial_guess = options->get(key+".initial_guess").as<string>();
             } else {
                 _initial_guess = "atom";
             }
            
            
             if ( options->exists(key+".externalfield_grid")) {
                 _do_externalfield=true;
                 _grid_name_ext=options->get(key + ".externalfield_grid").as<string>();
             }else{
                 _do_externalfield = false;
             }
            
            
	    // numerical integrations
	    _grid_name = options->get(key + ".integration_grid").as<string>();
            if ( options->exists(key+".integration_grid_small")) {
                 _use_small_grid=options->get(key+".integration_grid_small").as<bool>();
               }
               else{
                    _use_small_grid=true;
               }    
            _grid_name_small=Choosesmallgrid(_grid_name);
           
	    // exchange and correlation as in libXC
        
            _xc_functional_name = options->get(key + ".xc_functional").as<string> ();
        
	   
            _numofelectrons =0;
            
	   
            
            
            if ( options->exists(key+".convergence")) {
                
            if ( options->exists(key+".convergence.energy")) {
                  _Econverged=options->get(key+".convergence.energy").as<double>();
               }
               else{
                    _Econverged=1e-7;
               }    
            if ( options->exists(key+".convergence.error")) {
                  _error_converged=options->get(key+".convergence.error").as<double>();
               }
               else{
                    _error_converged=1e-7;
               }    
            
             if ( options->exists(key+".convergence.max_iterations")) {
                  _max_iter=options->get(key+".convergence.max_iterations").as<int>();
               }
               else{
                   _max_iter=100;
               }    
            
  
             if ( options->exists(key+".convergence.method")) {
               string method= options->get(key+".convergence.method").as<string>();
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
                   _useautomaticmixing=false;
                   _mixingparameter=options->get(key+".convergence.mixing").as<double>();
               }
               else{
                   _useautomaticmixing=true;
                    _mixingparameter=-10;
               }
            
             if ( options->exists(key+".convergence.levelshift")) {
                 
                   _levelshift=options->get(key+".convergence.levelshift").as<double>();
               }
               else{
                   _levelshift=0.0;
                    
               }
            
             if ( options->exists(key+".convergence.levelshift_end")) {
                 
                   _levelshiftend=options->get(key+".convergence.levelshift_end").as<double>();
               }
               else{
                   _levelshiftend=0.8;
                    
               }
            
            if ( options->exists(key+".convergence.DIIS_maxout")) {
                    _maxout=options->get(key+".convergence.DIIS_maxout").as<bool>();
               }
               else{
                     _maxout=false;
               }
             
              if ( options->exists(key+".convergence.DIIS_length")) {
                   _histlength=options->get(key+".convergence.DIIS_length").as<int>();
               }
               
               else{
                    _histlength=10;         
               }
               if ( options->exists(key+".convergence.DIIS_start")) {
                  _diis_start=options->get(key+".convergence.DIIS_start").as<double>();
               }
               else{
                    _diis_start=0.01;
               } 
             if ( options->exists(key+".convergence.ADIIS_start")) {
                  _adiis_start=options->get(key+".convergence.ADIIS_start").as<double>();
               }
               else{
                    _adiis_start=2;
               } 
        }
            else{
                _Econverged=1e-7;
                _error_converged=1e-7;
                 _maxout=false;
                 _diis_start=0.01;
                 _adiis_start=2;
                  _histlength=10;
                  _useautomaticmixing=true;
                  _mixingparameter=-10;
                   _usediis=true;
                   _max_iter=100;
                   _levelshift=0.25;
                   _levelshiftend=0.8;
            }
        if(!_usediis){
                _histlength=1;
                _maxout=false;                
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
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
            }
            #endif

            _atoms = _orbitals->QMAtoms();
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Molecule Coordinates [A] "  << flush;
            for(unsigned i=0;i<_atoms.size();i++){
                LOG(ctp::logDEBUG, *_pLog) << "\t\t "<< _atoms[i]->type<<" "<<_atoms[i]->x<<" "<<_atoms[i]->y<<" "<<_atoms[i]->z<<" "<<flush;
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

            ub::matrix<double> H0 = _dftAOkinetic.Matrix() + _dftAOESP.getNuclearpotential(); 
            
            NuclearRepulsion();
            if(_addexternalsites){
               H0+= _dftAOESP.getExternalpotential();
               //cout<<"analytic"<<_dftAOESP.getExternalpotential()<<endl;
               H0+= _dftAODipole_Potential.getExternalpotential();
               //H0+= _dftAOQuadrupole_Potential.getExternalpotential();
               
                double estat=ExternalRepulsion();
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " E_electrostatic "<<estat<<flush;
                E_nucnuc+=estat;
                
            }
            
             if(_do_externalfield){
                 //This was for testing purposes, now it works. :D
                
                 /*
            std::vector<const vec*> grid=_gridIntegration_ext.getGridpoints();
            for(std::vector<const vec*>::const_iterator gt=grid.begin();gt<grid.end();++gt){
                double value=0.0;
                for(std::vector<APolarSite*>::iterator it=_externalsites.begin();it<_externalsites.end();++it){
                    vec pos=(*it)->getPos();
                    double charge=(*it)->getQ00();
                    value+=charge/abs(pos*tools::conv::nm2bohr-*(*gt));
                }
                externalgrid.push_back(value);
            }
            
             std::vector<double> externalgrid_nuc;
             for(std::vector<::QMAtom*>::const_iterator at=_atoms.begin();at<_atoms.end();++at){
                double value=0.0;
                for(std::vector<APolarSite*>::iterator it=_externalsites.begin();it<_externalsites.end();++it){
                    vec pos=(*it)->getPos();
                    double charge=(*it)->getQ00();
                    value+=charge/abs(pos*tools::conv::nm2bohr-vec((*at)->x,(*at)->y,(*at)->z)*tools::conv::ang2bohr);
                }
                externalgrid_nuc.push_back(value);
            }
              */   
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Integrated external potential on grid "<< flush;
            //cout<<"grid"<<_gridIntegration_ext.IntegrateExternalPotential_Atomblock(&_dftbasis,externalgrid)<<endl;
            H0-=_gridIntegration_ext.IntegrateExternalPotential_Atomblock(&_dftbasis,_externalgrid);
            
            E_nucnuc+=ExternalGridRepulsion(_externalgrid_nuc);
            
            }
             
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Nuclear Repulsion Energy is " << E_nucnuc << flush;

            if(_with_ecp){
            H0+=_dftAOECP.Matrix();
            //cout<<_dftAOECP.Matrix()<<endl;
            cout<<endl;
            cout<< "WARNING ecps are up to numercis correct" <<endl;
            }
            
            // if we have a guess we do not need this. 
            if(!_with_guess){
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup Initial Guess using: "<<_initial_guess<< flush;
                 // this temp is necessary because eigenvalues_general returns MO^T and not MO
                if(_initial_guess=="independent"){
                _diis.SolveFockmatrix(MOEnergies,MOCoeff,H0);
                _dftAOdmat=_orbitals->DensityMatrixGroundState(MOCoeff);
                }
                else if (_initial_guess == "atom") {

                    _dftAOdmat = AtomicGuess(_orbitals);
                    MOEnergies=ub::zero_vector<double>(_dftAOdmat.size1());
                    //cout<<_dftAOdmat<<endl;
                    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<" Full atomic density Matrix gives N="<<std::setprecision(9)<<linalg_traceofProd(_dftAOdmat,_dftAOoverlap.Matrix())<<" electrons."<<flush;
                   
                    } else {
                        throw runtime_error("Initial guess method not known/implemented");
                    
                    }
                }else{
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Reading guess from orbitals object/file"<< flush;
                _dftbasis.ReorderMOs(MOCoeff, _orbitals->getQMpackage(), "xtp");
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converted DFT orbital coefficient order from " << _orbitals->getQMpackage() << " to xtp" << flush;
                _dftAOdmat=_orbitals->DensityMatrixGroundState(MOCoeff);
            }
           
            _orbitals->setQMpackage("xtp");
            
	    
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " STARTING SCF cycle" << flush;
            LOG(ctp::logDEBUG, *_pLog) << " --------------------------------------------------------------------------" << flush;
           
           double energyold=0;
           double diiserror=100;//is evolved in DIIs scheme
           Mixing Mixer(_useautomaticmixing,_mixingparameter,&_dftAOoverlap.Matrix(),  _pLog);
            for ( _this_iter=0; _this_iter<_max_iter; _this_iter++){
                LOG(ctp::logDEBUG, *_pLog)<< flush;
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Iteration "<< _this_iter+1 <<" of "<< _max_iter << flush;

                if(_with_RI){
                _ERIs.CalculateERIs(_dftAOdmat, _AuxAOcoulomb_inv);
                }
                else{
                    if(_4cmethod=="ram"){
                _ERIs.CalculateERIs_4c_small_molecule(_dftAOdmat);
                    }
               //     else if(_4cmethod=="direct"){
               // _ERIs.CalculateERIs_4c_large_molecule(_dftAOdmat);
                //    }
                }
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Electron repulsion matrix of dimension: " << _ERIs.getSize1() << " x " << _ERIs.getSize2()<<flush;
                double vxcenergy=0.0;
                if(_use_small_grid && diiserror>0.01){
                    _orbitals->AOVxc()=_gridIntegration_small.IntegrateVXC_Atomblock(_dftAOdmat,  &_dftbasis,_xc_functional_name);
                    vxcenergy=_gridIntegration_small.getTotEcontribution();
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled approximate DFT Vxc matrix "<<flush;           
                }
                else{
		_orbitals->AOVxc()=_gridIntegration.IntegrateVXC_Atomblock(_dftAOdmat,  &_dftbasis,_xc_functional_name);
                vxcenergy=_gridIntegration.getTotEcontribution();
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Vxc matrix "<<flush;              
                }
                
                ub::matrix<double> H=H0+_ERIs.getERIs()+_orbitals->AOVxc();
                
                double Eone=linalg_traceofProd(_dftAOdmat,H0);
                double Etwo=0.5*_ERIs.getERIsenergy()+vxcenergy;
                double totenergy=Eone+E_nucnuc+Etwo;
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Single particle energy "<<std::setprecision(12)<<Eone<<flush;
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Two particle energy "<<std::setprecision(12)<<Etwo<<flush;
                
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<std::setprecision(12)<< " Exc contribution "<<vxcenergy<<flush;
               
                
                
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total Energy "<<std::setprecision(12)<<totenergy<<flush;
                //cout<<_orbitals->AOVxc()<<endl;
                
       
                //this updates the density matrix as well
                diiserror=_diis.Evolve(_dftAOdmat,H,MOEnergies,MOCoeff,_this_iter,totenergy);
        
                
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " DIIs error "<<diiserror << flush;
                
                     ub::matrix<double> dmatin=_dftAOdmat;
                     _dftAOdmat=_orbitals->DensityMatrixGroundState(MOCoeff);
                     
                     ub::matrix<double>bla=ub::prod(_dftAOoverlap.Matrix(),ub::trans(MOCoeff));
                     ub::matrix<double>bla2=ub::prod(ub::trans(bla),_dftAOdmat);
                     ub::matrix<double>bla3=ub::prod(bla2,bla);
                     cout<<bla3<<endl;
                if (!(diiserror<_adiis_start && _usediis && _this_iter>2)){
                    _dftAOdmat=Mixer.MixDmat(dmatin,_dftAOdmat); 
                            //cout<<"mixing_beta"<<endl;
                        }
                    else{
                    Mixer.Updatemix(dmatin,_dftAOdmat);
                    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Using DIIs "<<flush;
                 }
                
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Updated Density Matrix "<<flush;
                
                
                if(tools::globals::verbose){
                for (int i=0;i<int(MOEnergies.size());i++){
                    if ( i <= _numofelectrons/2-1) {
                        LOG(ctp::logDEBUG, *_pLog) <<"\t\t" << i <<  " occ " << MOEnergies(i)  << flush;     
                    } else {
                        LOG(ctp::logDEBUG, *_pLog) <<"\t\t"<<   i <<  " vir " << MOEnergies(i)  << flush;
                        
                    }
                }
                }
                
                LOG(ctp::logDEBUG, *_pLog) << "\t\tGAP " << MOEnergies(_numofelectrons/2)-MOEnergies(_numofelectrons/2-1) << flush;
                
                 
                
                
              //  LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Solved general eigenproblem "<<flush;
                if (std::abs(totenergy-energyold)< _Econverged && diiserror<_error_converged){
                    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Total Energy has converged to "<<std::setprecision(9)<<std::abs(totenergy-energyold)<<"[Ha] after "<< _this_iter+1<<
                            " iterations. DIIS error is converged up to "<<_error_converged<<"[Ha]" <<flush;
                     LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Final Single Point Energy "<<std::setprecision(12)<<totenergy<<" Ha"<<flush;
                      LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " MO Energies  [Ha]"<<flush;
                      for (int i=0;i<int(MOEnergies.size());i++){
                    if ( i <= _numofelectrons/2-1) {
                        LOG(ctp::logDEBUG, *_pLog) <<"\t\t" << i <<  " occ "<<std::setprecision(12) << MOEnergies(i)  << flush;     
                    } else {
                        LOG(ctp::logDEBUG, *_pLog) <<"\t\t"<<   i <<  " vir "<<std::setprecision(12) << MOEnergies(i)  << flush;
                        
                    }
                }
                    break;
                }
                else{
                    energyold=totenergy;
                }
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<" Density Matrix gives N="<<std::setprecision(9)<<linalg_traceofProd(_dftAOdmat,_dftAOoverlap.Matrix())<<" electrons."<<flush;
                //LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Num of electrons "<< _gridIntegration.IntegrateDensity_Atomblock(_dftAOdmat, basis) << flush;
                
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
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftAOoverlap.Dimension() << flush;
            //cout<<"overlap"<<_dftAOoverlap.Matrix()<<endl;
	    // check DFT basis for linear dependence
            linalg_eigenvalues(_dftAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            
            //Not brilliant but we need S-1/2 for DIIS and I do not want to calculate it each time
            ub::matrix<double> _diagS = ub::zero_matrix<double>(_eigenvectors.size1(),_eigenvectors.size1() );
            for ( unsigned _i =0; _i < _eigenvalues.size() ; _i++){

                _diagS(_i,_i) = 1.0/sqrt(_eigenvalues[_i]);
            }
            ub::matrix<double> _temp = ub::prod( _diagS, ub::trans(_eigenvectors));
             _Sminusonehalf = ub::prod(_eigenvectors,_temp );
            
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Smallest eigenvalue of DFT Overlap matrix : " << _eigenvalues[0] << flush;
      }
            _dftAOkinetic.Initialize(_dftbasis.AOBasisSize());
            _dftAOkinetic.Fill(&_dftbasis);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT Kinetic energy matrix of dimension: " << _dftAOkinetic.Dimension() << flush;

            
            _dftAOESP.Initialize(_dftbasis.AOBasisSize());
            _dftAOESP.Fillnucpotential(&_dftbasis, _atoms,_with_ecp);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT nuclear potential matrix of dimension: " << _dftAOESP.Dimension() << flush;
            //_dftAOESP.Print("NUC");

            
            if (_addexternalsites){
                _dftAOESP.Fillextpotential(&_dftbasis, _externalsites);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external pointcharge potential matrix of dimension: " << _dftAOESP.Dimension() << flush;
                
                _dftAODipole_Potential.Fillextpotential(&_dftbasis, _externalsites);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external dipole potential matrix of dimension: " << _dftAODipole_Potential.Dimension() << flush;
                _dftAOQuadrupole_Potential.Fillextpotential(&_dftbasis, _externalsites);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT external quadrupole potential matrix of dimension: " << _dftAOQuadrupole_Potential.Dimension() << flush;
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " External sites\t Name \t Coordinates \t charge \t dipole \t quadrupole" << flush;
                for(unsigned i=0;i<_externalsites.size();i++){
                LOG(ctp::logDEBUG, *_pLog) << "\t\t "<< _externalsites[i]->getName()<<" | "<< _externalsites[i]->getPos().getX()
                                        <<" "<<_externalsites[i]->getPos().getY()<<" "<<_externalsites[i]->getPos().getZ()
                                        <<" | "<<_externalsites[i]->getQ00()<<" | "<<_externalsites[i]->getQ1().getX()
                                        <<" "<<_externalsites[i]->getQ1().getY()<<" "<<_externalsites[i]->getQ1().getZ()<<" | "
                                        <<_externalsites[i]->getQ2()[0]<<" "<<_externalsites[i]->getQ2()[1]<<" "<<_externalsites[i]->getQ2()[2]<<" "
                                                <<_externalsites[i]->getQ2()[3]<<" "<<_externalsites[i]->getQ2()[4]<<flush;
            }
                
            }
            //this will not remain here but be moved to qmape
            if(_do_externalfield){
                _gridIntegration_ext.GridSetup(_grid_name_ext,&_dftbasisset,_atoms,&_dftbasis);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name_ext 
                        << " for external field with "<<_gridIntegration_ext.getGridpoints().size()<<" points"<< flush;
            }
            
            if (_with_ecp) {
                _dftAOECP.Initialize(_dftbasis.AOBasisSize());
                _dftAOECP.Fill(&_dftbasis, vec(0,0,0), &_ecp);
                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled DFT ECP matrix of dimension: " << _dftAOECP.Dimension() << flush;
                //_dftAOECP.Print("ECP");
                
                _dftAOESP.getNuclearpotential() += _dftAOECP.Matrix();
                
            }
             
            
            _diis.Configure(_usediis, _histlength, _maxout, _diismethod, _adiis_start, _diis_start,_levelshift,_levelshiftend,false, _numofelectrons/2);
            _diis.setLogger(_pLog);
            _diis.setOverlap(&_dftAOoverlap.Matrix());
            _diis.setSqrtOverlap(&_Sminusonehalf);
	    // AUX AOoverlap
        
            if(_with_RI){
            { // this is just for info and not needed afterwards
            AOOverlap _auxAOoverlap;
            _auxAOoverlap.Initialize(_auxbasis.AOBasisSize());
            _auxAOoverlap.Fill(&_auxbasis);
            
            
            linalg_eigenvalues(_auxAOoverlap.Matrix(), _eigenvalues, _eigenvectors);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Smallest eigenvalue of AUX Overlap matrix : " << _eigenvalues[0] << flush;
             }
            
            
            
            
            
            // AUX AOcoulomb matrix
            //cout << " _auxAOcoulomb" <<endl;
            {
            AOCoulomb                           _auxAOcoulomb;
            _auxAOcoulomb.Initialize(_auxbasis.AOBasisSize());
            _auxAOcoulomb.Fill(&_auxbasis);
           
            //cout << _auxAOcoulomb._aomatrix<<endl;
            // _auxAOcoulomb.Print("COU");
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled AUX Coulomb matrix of dimension: " << _auxAOcoulomb.Dimension() << flush;
            _AuxAOcoulomb_inv=ub::zero_matrix<double>( _auxAOcoulomb.Dimension(), _auxAOcoulomb.Dimension()); 
            ub::matrix<double> _AuxAOcoulomb_inv2=ub::zero_matrix<double>( _auxAOcoulomb.Dimension(), _auxAOcoulomb.Dimension()); 
            int dimensions=linalg_invert_svd( _auxAOcoulomb.Matrix() , _AuxAOcoulomb_inv,1e8);
            
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Inverted AUX Coulomb matrix, removed "<<dimensions<<" functions from aux basis" << flush;
     
            }
            // prepare invariant part of electron repulsion integrals
      
            _ERIs.Initialize(_dftbasis, _auxbasis);
          
            
            
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup invariant parts of Electron Repulsion integrals " << flush;
            }
            else{
                if(_4cmethod=="ram"){
               _ERIs.Initialize_4c_small_molecule(_dftbasis); 
               LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated 4c integrals. " << flush;
                }
            }

      }
      
ub::matrix<double> DFTENGINE::AtomicGuess(Orbitals* _orbitals) {
            ub::matrix<double> guess = ub::zero_matrix<double>(_dftbasis.AOBasisSize());
            const std::vector<ctp::QMAtom*> atoms = _orbitals->QMAtoms();
            std::vector<ctp::QMAtom*> uniqueelements;
            std::vector<ctp::QMAtom*>::const_iterator at;
            std::vector<ctp::QMAtom*>::iterator st;
            std::vector< ub::matrix<double> > uniqueatom_guesses;
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Scanning molecule of size " << atoms.size() << " for unique elements" << flush;
            for (at = atoms.begin(); at < atoms.end(); ++at) {
                bool exists = false;
                if (uniqueelements.size() == 0) {
                    exists = false;
                } else {
                    for (st = uniqueelements.begin(); st < uniqueelements.end(); ++st) {
                        if ((*at)->type == (*st)->type) {
                            exists = true;
                            break;
                        }
                    }
                }
                if (!exists) {
                    uniqueelements.push_back((*at));
                }
            }
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " " << uniqueelements.size() << " unique elements found" << flush;
            Elements _elements;
            for (st = uniqueelements.begin(); st < uniqueelements.end(); ++st) {

                

                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculating atom density for " << (*st)->type << flush;
                bool with_ecp = _with_ecp;
                if ((*st)->type == "H" || (*st)->type == "He") {
                    with_ecp = false;
                }
                std::vector<ctp::QMAtom*> atom;
                atom.push_back(*st);

                AOBasis dftbasis;
                AOBasis ecp;
                NumericalIntegration gridIntegration;
                dftbasis.AOBasisFill(&_dftbasisset, atom);
                if (with_ecp) {
                    ecp.ECPFill(&_ecpbasisset, atom);
                }
                gridIntegration.GridSetup(_grid_name, &_dftbasisset, atom,&dftbasis);

                int numofelectrons = int(_elements.getNucCrg((*st)->type));
                int alpha_e=0;
                int beta_e=0;
                if (with_ecp) {
                    numofelectrons -= int(_ecpbasisset.getElement((*st)->type)->getNcore());
                }

               
                if ((numofelectrons%2)!=0){
                    alpha_e=numofelectrons/2+numofelectrons%2;
                    beta_e=numofelectrons/2;
                }
                else{
                    alpha_e=numofelectrons/2;
                    beta_e=alpha_e;
                }
                //cout<<"alpha "<<alpha_e<<endl;
                //cout<<"beta "<<beta_e<<endl;
              
             
                AOOverlap dftAOoverlap;
                AOKinetic dftAOkinetic;
                AOESP dftAOESP;
                AOECP dftAOECP;
                ub::vector<double> eigenvalues;
                ub::matrix<double> eigenvectors;
                ub::matrix<double> Sminusonehalf;
                ERIs ERIs_atom;


                // DFT AOOverlap matrix
                dftAOoverlap.Initialize(dftbasis.AOBasisSize());
                dftAOoverlap.Fill(&dftbasis);


                linalg_eigenvalues(dftAOoverlap.Matrix(), eigenvalues, eigenvectors);

                //Not brilliant but we need S-1/2 for DIIS and I do not want to calculate it each time
                ub::matrix<double> _diagS = ub::zero_matrix<double>(eigenvectors.size1(), eigenvectors.size1());
                for (unsigned _i = 0; _i < eigenvalues.size(); _i++) {

                    _diagS(_i, _i) = 1.0 / sqrt(eigenvalues[_i]);
                }
                ub::matrix<double> _temp = ub::prod(_diagS, ub::trans(eigenvectors));
                Sminusonehalf = ub::prod(eigenvectors, _temp);



                dftAOkinetic.Initialize(dftbasis.AOBasisSize());
                dftAOkinetic.Fill(&dftbasis);
                dftAOESP.Initialize(dftbasis.AOBasisSize());
                dftAOESP.Fillnucpotential(&dftbasis, atom, with_ecp);
                
                
                    
                  
                ERIs_atom.Initialize_4c_small_molecule(dftbasis);

                ub::vector<double>MOEnergies_alpha;
                ub::matrix<double>MOCoeff_alpha;
                ub::vector<double>MOEnergies_beta;
                ub::matrix<double>MOCoeff_beta;
                Diis diis_alpha;
                Diis diis_beta;
                
                Mixing Mix_alpha(false,0.7,&dftAOoverlap.Matrix(),_pLog);
                Mixing Mix_beta(false,0.7,&dftAOoverlap.Matrix(),_pLog);
                diis_alpha.Configure(true, 20, 0, "",  0.01,0.01,0.2,0.2,true,alpha_e );
                diis_alpha.setLogger(_pLog);
                diis_alpha.setOverlap(&dftAOoverlap.Matrix());
                diis_alpha.setSqrtOverlap(&Sminusonehalf);
                diis_beta.Configure(true, 20, 0, "",  0.01,0.01,0.2,0.2,true,beta_e );
                diis_beta.setLogger(_pLog);
                diis_beta.setOverlap(&dftAOoverlap.Matrix());
                diis_beta.setSqrtOverlap(&Sminusonehalf);
                /**** Construct initial density  ****/
               
                ub::matrix<double> H0 = dftAOkinetic.Matrix() + dftAOESP.getNuclearpotential();
                if (with_ecp) {
                    dftAOECP.Initialize(dftbasis.AOBasisSize());
                    dftAOECP.Fill(&dftbasis, vec(0,0,0), &ecp);
                    
                    H0 += dftAOECP.Matrix();
                    
                }
                
               
                diis_alpha.SolveFockmatrix(MOEnergies_alpha,MOCoeff_alpha,H0);
 
               MOEnergies_beta=MOEnergies_alpha;
               MOCoeff_beta=MOCoeff_alpha;
               cout<<MOEnergies_alpha<<endl;
            
               //ub::matrix<double>dftAOdmat_alpha = DensityMatrix_frac(MOCoeff_alpha,MOEnergies_alpha,alpha_e);
               //ub::matrix<double>dftAOdmat_beta = DensityMatrix_frac(MOCoeff_beta,MOEnergies_beta,beta_e);
               ub::matrix<double>dftAOdmat_alpha = DensityMatrix_unres(MOCoeff_alpha,alpha_e);
               ub::matrix<double>dftAOdmat_beta = DensityMatrix_unres(MOCoeff_beta,beta_e);

                    
              
               /* 
                double e_a=0.0;
                double e_b=0.0;
                 for (unsigned i = 0; i < MOEnergies_alpha.size(); i++) {
                       e_a+=occupation_alpha(i);
                   }
                   for (unsigned i = 0; i < MOEnergies_beta.size(); i++) {
                       e_b+=occupation_beta(i);
                   }
                cout <<"e_a "<<e_a<<endl;
                cout <<"e_b "<<e_b<<endl;
                 */ 
                    double energyold = 0;
                    int maxiter=201;
                    for (int this_iter = 0; this_iter < maxiter; this_iter++) {
                        //cout<<this_iter<<endl;
                        ERIs_atom.CalculateERIs_4c_small_molecule(dftAOdmat_alpha+dftAOdmat_beta);
                        ub::matrix<double> AOVxc_alpha = gridIntegration.IntegrateVXC_Atomblock(dftAOdmat_alpha, &dftbasis, _xc_functional_name);
                        //cout<<AOVxc_alpha<<endl;
                        double E_vxc_alpha= gridIntegration.getTotEcontribution();
                        ub::matrix<double> AOVxc_beta = gridIntegration.IntegrateVXC_Atomblock(dftAOdmat_beta, &dftbasis, _xc_functional_name);
                        double E_vxc_beta= gridIntegration.getTotEcontribution();
                        ub::matrix<double> H_alpha = H0 + ERIs_atom.getERIs() + AOVxc_alpha;
                        ub::matrix<double> H_beta = H0 + ERIs_atom.getERIs() + AOVxc_beta;
                        
                        
             
                        double E_one_alpha=linalg_traceofProd(dftAOdmat_alpha,H0);
                        double E_two_alpha=E_vxc_alpha+linalg_traceofProd(ERIs_atom.getERIs(),dftAOdmat_alpha);
                        double E_one_beta=linalg_traceofProd(dftAOdmat_beta,H0);
                        double E_two_beta=E_vxc_beta+linalg_traceofProd(ERIs_atom.getERIs(),dftAOdmat_beta);
                        double E_alpha=E_one_alpha+E_two_alpha;
                        double E_beta=E_one_beta+E_two_beta;
                   
                        
                        
                        double totenergy = E_alpha+E_beta;
                        //evolve alpha
                      double diiserror_alpha=diis_alpha.Evolve(dftAOdmat_alpha,H_alpha,MOEnergies_alpha,MOCoeff_alpha,this_iter,E_alpha);
                       
                            ub::matrix<double> dmatin=dftAOdmat_alpha;
                           dftAOdmat_alpha=DensityMatrix_unres(MOCoeff_alpha,alpha_e);
                           //dftAOdmat_alpha=DensityMatrix_frac(MOCoeff_alpha,MOEnergies_alpha,alpha_e);
                             if (!(diiserror_alpha<_adiis_start && _usediis && this_iter>2)){
                                  dftAOdmat_alpha=Mix_alpha.MixDmat(dmatin,dftAOdmat_alpha,false); 
                            //cout<<"mixing_alpha"<<endl;
                        }
                        else{
                         Mix_alpha.Updatemix(dmatin,dftAOdmat_alpha);
                        } 
                      //evolve beta
                      
                       double diiserror_beta=0.0;
                       if(beta_e>0){
                       diiserror_beta=diis_beta.Evolve(dftAOdmat_beta,H_beta,MOEnergies_beta,MOCoeff_beta,this_iter,E_beta);
                            ub::matrix<double> dmatin=dftAOdmat_beta;
                           dftAOdmat_beta=DensityMatrix_unres(MOCoeff_beta,beta_e);
                           //dftAOdmat_beta=DensityMatrix_frac(MOCoeff_beta,MOEnergies_beta,beta_e);
                        if (!(diiserror_beta<_adiis_start && _usediis && this_iter>2)){   
                            dftAOdmat_beta=Mix_beta.MixDmat(dmatin,dftAOdmat_beta,false); 
                            //cout<<"mixing_beta"<<endl;
                        }
                        else{
                            Mix_beta.Updatemix(dmatin,dftAOdmat_beta);
                        } 
                       }
                       
                  
                       cout<<MOEnergies_alpha<<endl;
                       cout<<MOEnergies_alpha(alpha_e)<<" "<<MOEnergies_alpha(alpha_e-1)<<endl;
                        
                         if(tools::globals::verbose){
                        LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()<<" Iter "<<this_iter<<" of "<<maxiter<<" Etot "<<totenergy<<" diise_a "<<diiserror_alpha<<" diise_b "<<diiserror_beta
                                <<" a_gap "<<MOEnergies_alpha(alpha_e)-MOEnergies_alpha(alpha_e-1)<<" b_gap "<<MOEnergies_beta(beta_e)-MOEnergies_beta(beta_e-1)<<flush;
                         }
                       bool converged=(std::abs(totenergy - energyold) < _Econverged && diiserror_alpha < _error_converged && diiserror_beta < _error_converged);
                        if (converged || this_iter==maxiter-1) {
                            
                            if(converged){
                                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << this_iter<<" iterations" << flush;
                            }
                            else{
                                LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Not converged after " << this_iter<<" iterations. Using unconverged density." << flush;
                            }
                            
                            
                            //ub::matrix<double> avdmat=AverageShells(dftAOdmat_alpha+dftAOdmat_beta,dftbasis);
                            
                            uniqueatom_guesses.push_back(dftAOdmat_alpha+dftAOdmat_beta);
                            //LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<" Atomic density Matrix for "<< (*st)->type<<" gives N="<<std::setprecision(9)<<linalg_traceofProd(avdmat,dftAOoverlap.Matrix())<<" electrons."<<flush;
                            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<" Atomic density Matrix for "<< (*st)->type<<" gives N="<<std::setprecision(9)<<linalg_traceofProd(dftAOdmat_alpha+dftAOdmat_beta,dftAOoverlap.Matrix())<<" electrons."<<flush;
                            break;
                  
                        }
                        else {
                            energyold = totenergy;
                        }


                    }
            }
            unsigned start=0;
            unsigned end=0;
            for (at = atoms.begin(); at < atoms.end(); ++at) {
                unsigned index=0;
                for (unsigned i=0; i<uniqueelements.size();i++) {
                    if((*at)->type==uniqueelements[i]->type){
                        index=i;
                        break;
                    }
                }
                
                Element* element = _dftbasisset.getElement((*at)->type);
                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                    end+=(*its)->getnumofFunc();
                }
                ub::project(guess, ub::range ( start, end ),ub::range (start, end )  )=uniqueatom_guesses[index];
                
                
                start=end;
            }
            
            return guess;
        }


      // PREPARATION 
      void DFTENGINE::Prepare( Orbitals* _orbitals ){

            // load and fill DFT basis set
            _dftbasisset.LoadBasisSet(_dftbasis_name);
            
            if(_with_guess){
                
                if (_orbitals->hasDFTbasis()){
                    if( _orbitals->getDFTbasis()!=_dftbasis_name){
                    throw runtime_error((boost::format("Basisset Name in guess orb file and in dftengine option file differ %1% vs %2%") %_orbitals->getDFTbasis() %_dftbasis_name).str() );
                    }
                }else{
                    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " WARNING: Orbital file has no basisset information,using it as a guess might work or not for calculation with " << _dftbasis_name << flush;
                }        
            }
            
            _orbitals->setDFTbasis( _dftbasis_name );    
	    _dftbasis.AOBasisFill( &_dftbasisset, _atoms);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

            if(_with_RI){
	    // load and fill AUX basis set
            _auxbasisset.LoadBasisSet(_auxbasis_name);
            //_orbitals->setDFTbasis( _dftbasis_name );
	    _auxbasis.AOBasisFill( &_auxbasisset, _atoms);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded AUX Basis Set " << _auxbasis_name << flush;
            }
            if(_with_ecp){
            // load ECP (element-wise information) from xml file
            _ecpbasisset.LoadPseudopotentialSet(_ecp_name);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Loaded ECP library " <<  _ecp_name << flush;

            // fill auxiliary ECP basis by going through all atoms
            _ecp.ECPFill(&_ecpbasisset, _atoms);
            LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Filled ECP Basis of size " << _ecp._aoshells.size() << flush;
            }
            
	    // setup numerical integration grid
            _gridIntegration.GridSetup(_grid_name,&_dftbasisset,_atoms,&_dftbasis);
            
            
            
	    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup numerical integration grid " << _grid_name << " for vxc functional " 
                    << _xc_functional_name<<" with " <<_gridIntegration.getGridpoints().size()<<" points"<< flush;
            
            if(_use_small_grid){
                 _gridIntegration_small.GridSetup(_grid_name_small,&_dftbasisset,_atoms,&_dftbasis);

	    LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Setup small numerical integration grid " << _grid_name_small << " for vxc functional " 
                    << _xc_functional_name<<" with " <<_gridIntegration_small.getGridpoints().size()<<" points"<< flush;
            }
            
            
	
           Elements _elements; 
            //set number of electrons and such
           _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
           
           for (unsigned i=0;i<_atoms.size();i++){
               _numofelectrons+=_elements.getNucCrg(_atoms[i]->type);
           }

            // if ECP
            if (_with_ecp) {
                for (unsigned i = 0; i < _atoms.size(); i++) {
                    if (_atoms[i]->type=="H" || _atoms[i]->type=="He"){
                        continue;
                    }
                    else{
                    _numofelectrons-= _ecpbasisset.getElement(_atoms[i]->type)->getNcore() ;
                    }
                }
            }
           // here number of electrons is actually the total number, everywhere else in votca it is just alpha_electrons
           
           LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Total number of electrons: " << _numofelectrons << flush;
           
           if(_with_guess){
               if(_orbitals->getNumberOfElectrons()!=_numofelectrons/2){
               throw runtime_error((boost::format("Number of electron in guess orb file %1% and in dftengine differ %2%.")%_orbitals->getNumberOfElectrons() %(_numofelectrons/2) ).str());
               }
               if(_orbitals->getNumberOfLevels()!=_dftbasis.AOBasisSize()){
               throw runtime_error((boost::format("Number of levels in guess orb file %1% and in dftengine differ %2%.")%_orbitals->getNumberOfLevels() %_dftbasis.AOBasisSize()).str());
               }
            }
           else{
           _orbitals->setNumberOfElectrons(_numofelectrons/2);
           _orbitals->setNumberOfLevels(_numofelectrons/2,_dftbasis.AOBasisSize()-_numofelectrons/2);
           }
           

            // _orbitals->setBasisSetSize(_dftbasis.AOBasisSize());
            
      }
      
      
      
      
      ub::matrix<double> DFTENGINE::DensityMatrix_unres( const ub::matrix<double>& MOs, int numofelec ) { 
          if(numofelec==0){
              return ub::zero_matrix<double>(MOs.size1());
          }
         
         ub::matrix<double> _dmatGS = ub::zero_matrix<double>(MOs.size1());
        #pragma omp parallel for
        for ( unsigned _i=0; _i < MOs.size1(); _i++ ){
            for ( unsigned _j=0; _j < MOs.size1(); _j++ ){
                for ( int _level=0; _level < numofelec; _level++ ){
                 
                    _dmatGS(_i,_j) +=  MOs( _level , _i ) * MOs( _level , _j );
                 
                }
            }
         }
     //}    
     // return     
     return _dmatGS;  
 }
      
      
   ub::matrix<double> DFTENGINE::DensityMatrix_frac( const ub::matrix<double>& MOs,const ub::vector<double>& MOEnergies, int numofelec ) { 
          if(numofelec==0){
              return ub::zero_matrix<double>(MOs.size1());
          }
          
         ub::vector<double>occupation=ub::zero_vector<double>(MOEnergies.size());
         
         double buffer=0.0001;
         double homo_energy=MOEnergies(numofelec-1);
         std::vector<unsigned> degeneracies;
         
         for ( unsigned _level=0; _level < occupation.size(); _level++ ){
            if (MOEnergies(_level)<(homo_energy-buffer)){
                   occupation(_level)=1.0;
                   numofelec--;
               }
               else if(std::abs(MOEnergies(_level)-homo_energy)<buffer){
                   degeneracies.push_back(_level);
               }
               else if(MOEnergies(_level)>(homo_energy+buffer)){
                   occupation(_level)=0.0;
               }
           } 
         double deg_occupation=double(numofelec)/double(degeneracies.size());
         for ( unsigned _level=0; _level < degeneracies.size(); _level++ ){
             occupation(degeneracies[_level])=deg_occupation;         
         }
        cout<<occupation<<endl;
        ub::matrix<double> _dmatGS = ub::zero_matrix<double>(MOs.size1());
        #pragma omp parallel for
        for ( unsigned _i=0; _i < MOs.size1(); _i++ ){
            for ( unsigned _j=0; _j < MOs.size1(); _j++ ){
                for ( unsigned _level=0; _level < occupation.size(); _level++ ){
                 
                    _dmatGS(_i,_j) +=occupation(_level)*  MOs( _level , _i ) * MOs( _level , _j );
                 
                }
            }
         }
     //}    
     // return     
     return _dmatGS;  
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
      
      
      double DFTENGINE::ExternalRepulsion(){
          Elements element;
          double E_ext=0.0;
          
          if(_externalsites.size()==0){
              return 0;
          }
          
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
              for(std::vector<ctp::APolarSite*>::iterator ext=_externalsites.begin();ext<_externalsites.end();++ext){
                  vec r2=(*ext)->getPos()*tools::conv::nm2bohr;
                  double charge2=(*ext)->getQ00();
                  E_ext+=charge1*charge2/(abs(r1-r2));
                  if ((*ext)->getRank()>0 || (*ext)->IsPolarizable()){
                      vec dipole=((*ext)->getU1()+(*ext)->getQ1())*tools::conv::nm2bohr;
                      E_ext-=charge1*dipole*(r2-r1)/pow(abs(r2-r1),3);
                      
                      
                  }
                  if((*ext)->getRank()>1){
                      cout <<endl;
                      cout<<"WARNING: external multipoles higher than dipoles are not yet taken into account!"<<endl;
                  }
                  
              }
          }
          
          return E_ext;
      }
      
      
      double DFTENGINE::ExternalGridRepulsion(std::vector<double> externalpotential_nuc){
          Elements element;
          double E_ext=0.0;
          
          if(!_do_externalfield){
              return 0;
          }

          for(unsigned i=0;i<_atoms.size();i++){
              double charge=0.0;
               if(_with_ecp){
                  charge=(element.getNucCrgECP(_atoms[i]->type));
              }
              else{
                  charge=(element.getNucCrg(_atoms[i]->type));
              }
                  E_ext+=charge*externalpotential_nuc[i];
              }
          

          return E_ext;
      }
      
      
    
      
      string DFTENGINE::Choosesmallgrid(string largegrid){
          string smallgrid;
          
          if (largegrid=="xfine"){
            smallgrid="fine";  
          }
          else if(largegrid=="fine"){
            smallgrid="medium";  
          }
          else if(largegrid=="medium"){
            _use_small_grid=false;
            smallgrid="medium";  
          }
          else if(largegrid=="coarse"){
             _use_small_grid=false;
            smallgrid="coarse";  
          }
          else if(largegrid=="xcoarse"){
            _use_small_grid=false;
            smallgrid="xcoarse";  
          }
          else{
               throw runtime_error("Grid name for Vxc integration not known.");
          }
          
          return smallgrid;
      }
      
      
      //average atom densities matrices, for SP and other combined shells average each subshell separately. Does not really work yet
      ub::matrix<double> DFTENGINE::AverageShells(const ub::matrix<double>& dmat, AOBasis& dftbasis){
          ub::matrix<double> avdmat=ub::zero_matrix<double>(dmat.size1());
          AOBasis::AOShellIterator it;
          AOBasis::AOShellIterator jt;
          int start1=0.0;
          for(it=dftbasis.firstShell();it<dftbasis.lastShell();++it){
            AOShell* shell1=dftbasis.getShell(it);
            int end1=shell1->getNumFunc()+start1;
            std::vector<int> starts1;
            std::vector<int> ends1;
            //cout<<shell1->getType()<<" Start1 "<<start1<<" End1 "<<end1<<endl;
            //check if shell is combined
            if(shell1->getLmax()!=shell1->getLmin()){
                std::vector<int> temp1=NumFuncSubShell(shell1->getType());
                int numfunc1=start1;
                for(unsigned i=0;i<temp1.size();i++){
                    
                    starts1.push_back(numfunc1);
                    numfunc1+=temp1[i];
                    ends1.push_back(numfunc1);
                }
            }
            else{
                starts1.push_back(start1);
                ends1.push_back(end1);
            }
            start1=end1;
            int start2=0;
              for(jt=dftbasis.firstShell();jt<dftbasis.lastShell();++jt){
                AOShell* shell2=dftbasis.getShell(jt);
               
                int end2=shell2->getNumFunc()+start2;
                //cout<<shell2->getType()<<" Start2 "<<start2<<" End2 "<<end2<<endl;
                std::vector<int> starts2;
                std::vector<int> ends2;
                if(shell2->getLmax()!=shell2->getLmin()){
                    std::vector<int> temp2=NumFuncSubShell(shell2->getType());
                    int numfunc2=start2;
                    for(unsigned i=0;i<temp2.size();i++){

                        starts2.push_back(numfunc2);
                        numfunc2+=temp2[i];
                        ends2.push_back(numfunc2);
                    }
                }
                else{
                    starts2.push_back(start2);
                    ends2.push_back(end2);
                }
                start2=end2;
                for(unsigned k=0;k<starts1.size();k++){
                    int s1=starts1[k];
                    int e1=ends1[k];
                    for(unsigned l=0;l<starts2.size();l++){
                    int s2=starts2[l];
                    int e2=ends2[l];

                    double avg=0.0;
                    int count=0.0;

                    for(int i=s1;i<e1;++i){
                       for(int j=s2;j<e2;++j){
                               avg+=dmat(i,j);
                               count++;   
                          } 
                    }
                    avg=avg/double(count);
                    for(int i=s1;i<e1;++i){
                       for(int j=s2;j<e2;++j){
                               avdmat(i,j)=avg;   
                          } 
                      }
                    }
                }
            }     
        }  
        return avdmat;  
        }
    }}
