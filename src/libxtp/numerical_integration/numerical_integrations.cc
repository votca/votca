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
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
//for libxc
#include <votca/xtp/votca_config.h>

#include <votca/xtp/numerical_integrations.h>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/radial_euler_maclaurin_rule.h>
#include <votca/xtp/sphere_lebedev_rule.h>
#include <votca/xtp/aoshell.h>
#include <votca/tools/constants.h>

#ifdef LIBXC
#include <xc.h>
#endif
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/exchange_correlation.h>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/xtp/vxc_functionals.h>
#include <iterator>
#include <string>




namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        double NumericalIntegration::getExactExchange(const string _functional){
#ifdef LIBXC            
        
            double exactexchange=0.0;
            Vxc_Functionals map;
            std::vector<string> strs;
            
            boost::split(strs, _functional, boost::is_any_of(" "));
            if (strs.size()>2 ) {
                throw std::runtime_error("Too many functional names");
            }
            else if (strs.size()<1 ) {
                throw std::runtime_error("Specify at least one funcitonal");
            }
            
            for (unsigned i=0;i<strs.size();i++){
               
                int func_id = map.getID(strs[i]); 
                if (func_id<0){
                    exactexchange=0.0;
                    break;
                }
                xc_func_type func;
                if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", func_id);
                    exit(1);
                }
                if (exactexchange>0 && func.cam_alpha>0){
                    throw std::runtime_error("You have specified two functionals with exact exchange");
                }
                exactexchange+=func.cam_alpha;
            
                
            
            }
            return exactexchange;
            
#else
            return 0.0;
#endif
            
        }
        
        ub::matrix<double> NumericalIntegration::IntegrateExternalPotential_Atomblock(const std::vector<double>& Potentialvalues){
            if(_significant_atoms.size()<1){
                throw runtime_error("NumericalIntegration::IntegrateExternalPotential_Atomblock:significant atoms not found yet.");
            }
            ub::matrix<double> ExternalMat = ub::zero_matrix<double>(_basis->AOBasisSize(), _basis->AOBasisSize());
            
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            // separate storage for each thread
            std::vector< ub::matrix<double> > expot_thread;
            
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){     
                expot_thread.push_back( ub::zero_matrix<double>(_basis->AOBasisSize(), _basis->AOBasisSize()) );              
            }

            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;
                
           
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )

                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, _basis->AOBasisSize()); // TRY MORE USEFUL DATA
                          
                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                  
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                    
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){
                      
                         //   boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                            AOBasis::AOShellIterator _row = _atomshells[rowatom][ishell];
                            // for density, fill sub-part of AOatgrid
                            //ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc(), 0, 1);
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                            // (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_x, _grid[i][j].grid_y, _grid[i][j].grid_z);

                        
                            (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);
                        

                        }  // shell in atom
                
                
                    } // row shells 

                
                   // yeah storin potential values in a vector is weird but I did not want to cram it into gridpoint, because that will blow the structure more than necessary
                    ub::matrix<double> _addExt = 0.5*_grid[i][j].grid_weight * AOgrid* Potentialvalues[i*_grid[i].size()+j];

                    // combine/sum atom-block wise, only trigonal part, symmetrize later
                    // for each significant atom for this grid point
                    // parallelization only accesses atomblock information (_addXC, AOgrid -> XCmatblock), so no trouble with shared memory access )
                    // #pragma omp parallel for
                    for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {
                        
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                    
                        ub::matrix_range< ub::matrix<double> > _rowExt = ub::subrange( _addExt, 0 , 1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);    

                        for (unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size(); sigcol++) {
                            int colatom = _significant_atoms[i][j][sigcol];
                            // if (colatom > rowatom) break;

                            ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange( AOgrid, 0,1,  _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            
                            // update block reference of XCMAT
                            ub::matrix_range<ub::matrix<double> > _expotmatblock = ub::subrange( expot_thread[i_thread],_startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom] );
                            //_XCmatblock += ub::prod( _rowXC, ub::trans(_AOcol)  );
                            _expotmatblock += ub::prod( ub::trans(_rowExt), _AOcol  );

                            // update the other block
  
                        } // significant col
                    } // significant row 

                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid


            // sum thread matrices
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                #pragma omp parallel for
                for (unsigned _i = 0; _i < ExternalMat.size1(); _i++) {
                    //for (int _j = 0; _j <= _i; _j++) {
                        for (unsigned _j = 0; _j <ExternalMat.size2(); _j++) {
                    ExternalMat( _i, _j ) += expot_thread[i_thread](_i, _j);
                    }
                }
            }
            
            ExternalMat += ub::trans(ExternalMat);

            return ExternalMat;

        }
        
        
        
        ub::matrix<double> NumericalIntegration::IntegrateVXC_Atomblock(const ub::matrix<double>& _density_matrix,const string _functional){
            EXC = 0;
            if(_significant_atoms.size()<1){
                throw runtime_error("NumericalIntegration::IntegrateVXC_Atomblock:significant atoms not found yet.");
            }
            // TODO: switch XC functionals implementation from LIBXC to base own calculation
            ExchangeCorrelation _xc;
            Vxc_Functionals map;
            std::vector<string> strs;           
            boost::split(strs, _functional, boost::is_any_of(" "));
            int xfunc_id = 0;
            
#ifdef LIBXC
            bool _use_votca = false;
            bool _use_separate = false;
            int cfunc_id = 0;

            if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
                if (xfunc_id < 0) _use_votca = true;
            }

            else if (strs.size() == 2) {
                cfunc_id = map.getID(strs[0]);
                xfunc_id = map.getID(strs[1]);
                _use_separate = true;
            }
            else {
                cout<<"LIBXC "<<strs.size()<<endl;
                throw std::runtime_error("With LIBXC. Please specify one combined or an exchange and a correlation functionals");

            }
            xc_func_type xfunc; // handle for exchange functional
            xc_func_type cfunc; // handle for correlation functional
            if (!_use_votca){
            if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
                fprintf(stderr, "Functional '%d' not found\n", xfunc_id);
                exit(1);
            }
            
            xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
            if (xfunc.info->kind!=2 && !_use_separate){
                throw std::runtime_error("Your functional misses either correlation or exchange, please specify another functional, separated by whitespace");
            }
            
            if (_use_separate) {
                if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
                    fprintf(stderr, "Functional '%d' not found\n", cfunc_id);
                    exit(1);
                }
                xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED);
                xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED);
                if ((xfunc.info->kind+cfunc.info->kind)!=1){
                    throw std::runtime_error("Your functionals are not one exchange and one correlation");
                }
            }
            }
#else
         if (strs.size() == 1) {
                xfunc_id = map.getID(strs[0]);
            }   
         else {
                throw std::runtime_error("Running without LIBXC, Please specify one combined or an exchange and a correlation functionals");
         }
#endif
            
            //split dmat into atomsize protions so that access is faster later on
             #pragma omp parallel for
            for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
                for (unsigned colatom=0;colatom<=rowatom;colatom++){
            
            dmat_vector[rowatom][colatom] = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
            }
        }
            
           
            
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            // separate storage for each thread
            //same as for dmat for vxc mat
            std::vector< ub::matrix<double> > XCMAT_thread;
            std::vector<double> EXC_thread;
             for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                  EXC_thread.push_back(0.0);
                  
            }
            #pragma omp parallel for
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                
                for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
                for (unsigned colatom=0;colatom<_grid.size();colatom++){
                    xcmat_vector_thread[i_thread][rowatom][colatom]= ub::zero_matrix<double>(_blocksize[rowatom], _blocksize[colatom]);
                }
                }
            }
            #pragma omp parallel for
            for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
                for (unsigned colatom=0;colatom<_grid.size();colatom++){
                    xcmat_vector[rowatom][colatom]= ub::zero_matrix<double>(_blocksize[rowatom], _blocksize[colatom]);
                    }
                }
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;
                
               
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {

                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, _basis->AOBasisSize()); // TRY MORE USEFUL DATA
		    // get value of density gradient at each gridpoint
                 
                   ub::matrix<double> gradAOgrid = ub::zero_matrix<double>(3, _basis->AOBasisSize()); // for Gradients of AOs
                    
                   ub::matrix<double>  rho_mat = ub::zero_matrix<double>(1,1);
                   //ub::matrix<double> grad_rho = ub::zero_matrix<double>(3,1);
                    
                   ub::matrix<double> grad_rho = ub::zero_matrix<double>(1,3);
		    // evaluate AO Functions for all shells, NOW BLOCKWISE
                   
                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                       
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                    
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){
                      
                        
                            AOBasis::AOShellIterator _row = _atomshells[rowatom][ishell];
                            // for density, fill sub-part of AOatgrid
                           
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                         

                            // gradient of density
                        
                            ub::matrix_range< ub::matrix<double> > _gradAO = ub::subrange(gradAOgrid, 0, 3, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                           
                            (*_row)->EvalAOspace(_AOgridsub, _gradAO , _grid[i][j].grid_pos);
                           
                      
                        }  // shell in atom
                
                        ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);
                        ub::matrix<double> _tempgrad = ub::zero_matrix<double>(3,_blocksize[rowatom]);
                        
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                        // for each atom
                        // for all significant atoms of triangular matrix
                        for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                          
                            int colatom = _significant_atoms[i][j][sigcol];
                            if ( colatom > rowatom ) break;
                            
                            // get the already calculated AO values

                            ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                            ub::matrix_range< ub::matrix<double> > _gradAOgridcol = ub::subrange(gradAOgrid, 0, 3, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);

                          
                            const ub::matrix<double> & DMAT_here = dmat_vector[rowatom][colatom];
                            
                             if ( colatom == rowatom ){
                                _temp     += 0.5 * ub::prod( _AOgridcol, DMAT_here);
                                _tempgrad += 0.5 * ub::prod( _gradAOgridcol, DMAT_here);
                            } else {
                                
                                _temp     += ub::prod(  _AOgridcol, DMAT_here);
                                _tempgrad += ub::prod( _gradAOgridcol, DMAT_here);
                            }

                        } //col shells
                                             
                       
                        ub::matrix_range< ub::matrix<double> > _gradAOgridrow = ub::subrange(gradAOgrid, 0,3, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
                   

                        rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                        grad_rho += ub::prod(_temp, ub::trans(_gradAOgridrow)) +  ub::prod(_AOgridrow,ub::trans(_tempgrad)) ;

                    } // row shells 

                    double rho      = 2.0 * rho_mat(0,0);
              
                    
		    if ( rho < 1.e-15 ) continue; // skip the rest, if density is very small
                    grad_rho = 2.0 * grad_rho;
                                     
                    // get XC for this density_at_grid
                    double f_xc;      // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho;   // v_xc_rho(r) = df/drho
                    double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
                   
 #ifdef LIBXC                   
                    if (_use_votca) {
#endif                  
                        cout<<"Warning: VOTCA_PBE does give correct Vxc but incorrect E_xc"<<endl;
                        _xc.getXC(xfunc_id, rho, grad_rho(0, 0), grad_rho(0, 1), grad_rho(0, 2), f_xc, df_drho, df_dsigma);
#ifdef LIBXC
                    }                        // evaluate via LIBXC, if compiled, otherwise, go via own implementation

                    else {
                     

                        double sigma = ub::prod(grad_rho, ub::trans(grad_rho))(0, 0);

                        double exc[1];
                        double vsigma[1]; // libxc 
                        double vrho[1]; // libxc df/drho
                        switch (xfunc.info->family) {
                            case XC_FAMILY_LDA:
                                xc_lda_exc_vxc(&xfunc, 1, &rho, exc, vrho);
                                break;
                            case XC_FAMILY_GGA:
                            case XC_FAMILY_HYB_GGA:
                                xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                break;
                        }
                        f_xc = exc[0];
                        df_drho = vrho[0];
                        df_dsigma = vsigma[0];
                        if (_use_separate) {
                            // via libxc correlation part only
                            switch (cfunc.info->family) {
                                case XC_FAMILY_LDA:
                                    xc_lda_exc_vxc(&cfunc, 1, &rho, exc, vrho);
                                    break;
                                case XC_FAMILY_GGA:
                                case XC_FAMILY_HYB_GGA:
                                    xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, exc, vrho, vsigma);
                                    break;
                            }

                            f_xc += exc[0];
                            df_drho += vrho[0];
                            df_dsigma += vsigma[0];
                        }
                    }
#endif
                  
                    ub::matrix<double> _addXC = _grid[i][j].grid_weight * df_drho * AOgrid *0.5;

                    _addXC+=  2.0*df_dsigma * _grid[i][j].grid_weight * ub::prod(grad_rho,gradAOgrid);

                    // Exchange correlation energy
                    EXC_thread[i_thread] += _grid[i][j].grid_weight * rho * f_xc;
  
                    // combine/sum atom-block wise, only trigonal part, symmetrize later
                    // for each significant atom for this grid point
                    // parallelization only accesses atomblock information (_addXC, AOgrid -> XCmatblock), so no trouble with shared memory access )
                    // #pragma omp parallel for
                    for (unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size(); sigrow++) {
                        
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                    
                        const ub::matrix_range< ub::matrix<double> > _rowXC = ub::subrange( _addXC, 0 , 1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);    

                     
                        std::vector< ub::matrix<double> >& _XCmatblock = xcmat_vector_thread[i_thread][rowatom];
                        for (unsigned sigcol = 0; sigcol <_significant_atoms[i][j].size(); sigcol++) {
                            int colatom = _significant_atoms[i][j][sigcol];
                            
                            const ub::matrix_range< ub::matrix<double> > _AOcol = ub::subrange( AOgrid, 0,1,  _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);                         
                            _XCmatblock[colatom]+= ub::prod( ub::trans(_rowXC), _AOcol  );
  
                        } // significant col
                    } // significant row 

                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid


            // sum thread matrices
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                EXC += EXC_thread[i_thread];
                
            }
            
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                #pragma omp parallel for
                for (unsigned _i = 0; _i < xcmat_vector.size(); _i++) {
                    for (unsigned _j = 0; _j < xcmat_vector[_i].size(); _j++) {
                   
                      
                    xcmat_vector[_i][_j] += xcmat_vector_thread[i_thread][_i][_j];
                    }
                }
            }
             
             ub::matrix<double> XCMAT = ub::zero_matrix<double>(_basis->AOBasisSize(), _basis->AOBasisSize());
             
             #pragma omp parallel for
             for (unsigned rowatom=0;rowatom<xcmat_vector.size();rowatom++){
                for (unsigned colatom=0;colatom<xcmat_vector[rowatom].size();colatom++){
                    
            ub::subrange( XCMAT, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom], _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom])
                    =xcmat_vector[rowatom][colatom];
            
            }
        }
         
            XCMAT+=ub::trans(XCMAT);   
  
            return XCMAT;
        }
        

            
        double NumericalIntegration::IntegratePotential(const vec& rvector){
            
            double result = 0.0;
            
           if(density_set){
                for (unsigned i = 0; i < _grid.size(); i++) {
                for (unsigned j = 0; j < _grid[i].size(); j++) {
                    double dist=abs((_grid[i][j].grid_pos-rvector));
                    result -= _grid[i][j].grid_weight * _grid[i][j].grid_density/dist;
                    }
                }
            } 
           else{
               throw std::runtime_error("Density not calculated");
           }
            
            return result;   
        }
               
        
        void NumericalIntegration::FindsignificantAtoms(){
            
            int _atomindex = 0;
            int _Idx       = 0;
            int _size      = 0;
            
            for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                 
                if ( (*_row)->getIndex() == _atomindex ){
                    
                    _singleatom.push_back(_row);
                    _size += (*_row)->getNumFunc();
                       
                } else {
                    
                    // append _singleatom to _atomshells
                    _atomshells.push_back(_singleatom);
                    _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
                    // reset _singleatom
                    _singleatom.clear();
                    _size = (*_row)->getNumFunc();
                    _Idx       = (*_row)->getStartIndex();
                    _singleatom.push_back(_row);
                    _atomindex = (*_row)->getIndex();
                    
                }   
            }
            
            _atomshells.push_back(_singleatom);
            _startIdx.push_back( _Idx );
                    _blocksize.push_back(_size);
          
           
            // setup a list of min decay constants per atom
            // for every shell
            _atomindex = 0;
            double _decaymin = 1e7;
            vector< double > _minimal_decay;
            vector < vec > _positions;
            vec _localpos = (*_basis->firstShell())->getPos();
            for ( AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++   ) {
                               
                 if ( (*_row)->getIndex() == _atomindex ){
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         const AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->getDecay();
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives 
                     
                 } else {  // if shell belongs to the actual atom
                     // add to mininal_decay vector
                     _minimal_decay.push_back(_decaymin);
                     _positions.push_back( _localpos );
                     // reset counters
                     _decaymin = 1e7;
                     _localpos = (*_row)->getPos();

                     _atomindex++;
                     
                     // check all decay constants in this shell
                     for (AOShell::GaussianIterator itg = (*_row)->firstGaussian(); itg != (*_row)->lastGaussian(); itg++) {
                         const AOGaussianPrimitive* gaussian = *itg;
                         double _decay = gaussian->getDecay();
                         if (_decay < _decaymin) {
                             _decaymin = _decay;
                         } // decay min check
                     
                     } // Gaussian Primitives                                       
                 }
            } // all shells
                 
            // push final atom
            _minimal_decay.push_back(_decaymin);
            _positions.push_back( _localpos );
            
                          
             // for each gridpoint, check the value of exp(-a*(r-R)^2) < 1e-10
             //                             = alpha*(r-R)^2 >~ 20.7
            
            // each atomic grid
            for (unsigned i = 0; i < _grid.size(); i++) {
            
                vector< vector<int> > _significant_atoms_atomgrid;
                
                // each point of the atomic grid
                for (unsigned j = 0; j < _grid[i].size(); j++) {

                    vector<int> _significant_atoms_gridpoint;
                    const vec& grid=_grid[i][j].grid_pos;
                   
                    
                    // check all atoms
                    for ( unsigned iatom = 0 ; iatom < _minimal_decay.size(); iatom++){

                        vec dist = grid - _positions[iatom];
                        double distsq = dist*dist ;
                        
                        // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (_minimal_decay[iatom] * distsq) < 20.7 ){
                            _significant_atoms_gridpoint.push_back(iatom);
                        }
                        
                    } // check all atoms

                    _significant_atoms_atomgrid.push_back(  _significant_atoms_gridpoint );
                   
                } // all points of this atom grid
                
                _significant_atoms.push_back(_significant_atoms_atomgrid);
               
            } // atomic grids
              
       
            
            int total_grid =0;
            int significant_grid = 0;
            for ( unsigned i = 0; i < _significant_atoms.size(); i++ ){
                
                total_grid += _grid[i].size(); 
                
                for ( unsigned j = 0; j < _significant_atoms[i].size(); j++ ){
                    
                    int gridpointsize = _significant_atoms[i][j].size();
                    significant_grid += gridpointsize*(gridpointsize+1);
                    
                } 
            }
            int natoms = _grid.size();
          
            total_grid = total_grid * ( natoms*(natoms+1) ) / 2;
            
         for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
            std::vector< ub::matrix<double> > rowmatrix; 
                      for (unsigned colatom=0;colatom<=rowatom;colatom++){
                         rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[colatom],_blocksize[rowatom]));
                 }
            dmat_vector.push_back(rowmatrix);
         }
            
            
            
             // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

             
               for(int i=0;i<nthreads;i++){
               
            std::vector< std::vector< ub::matrix<double> > > matrix; 
              for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
            std::vector< ub::matrix<double> > rowmatrix; 
                      for (unsigned colatom=0;colatom<_grid.size();colatom++){
                          rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[rowatom],_blocksize[colatom])); 
                 }
           matrix.push_back(rowmatrix);
         } 
            xcmat_vector_thread.push_back(matrix);
               }
               
    for (unsigned rowatom=0;rowatom<_grid.size();rowatom++){
        std::vector< ub::matrix<double> > rowmatrix; 
          for (unsigned colatom=0;colatom<_grid.size();colatom++){
              rowmatrix.push_back(ub::zero_matrix<double>(_blocksize[colatom],_blocksize[rowatom]));
             }
       xcmat_vector.push_back(rowmatrix);
     } 
        return;
        }
        

        double NumericalIntegration::IntegrateDensity_Atomblock(const ub::matrix<double>& _density_matrix){   
            if(_significant_atoms.size()<1){
                throw runtime_error("NumericalIntegration::IntegrateDensity_Atomblock:significant atoms not found yet.");
            }
            double result=0.0;

             
          
            
  
            // parallelization: distribute over threads inside one atom
            int nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif

            std::vector<double> Density_thread;
            for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){ 
                Density_thread.push_back(0.0);
            }           
            
            // for every atom
            for (unsigned i = 0; i < _grid.size(); i++) {
	      // for each point in atom grid
                
                // number of points in this atomgrid
                int atom_points = _grid[i].size();
                // divide among threads
                int atom_points_per_thread = atom_points/nthreads;
                std::vector<int> _thread_start;
                std::vector<int> _thread_stop;
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                    _thread_start.push_back( i_thread * atom_points_per_thread );
                    _thread_stop.push_back( (i_thread + 1) * atom_points_per_thread );
                }
                // final stop must be size
                _thread_stop[nthreads-1] = atom_points;

         
                
                #pragma omp parallel for
                for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                for (int j = _thread_start[i_thread]; j < _thread_stop[i_thread]; j++) {
                   //boost::timer::cpu_times t0 = cpu_t.elapsed();

                    // get value of orbitals at each gridpoint (vector as 1D boost matrix object -> prod )
                    //ub::matrix<double> AOgrid = ub::zero_matrix<double>(basis->AOBasisSize(), 1);

                   ub::matrix<double> AOgrid = ub::zero_matrix<double>(1, _basis->AOBasisSize()); // TRY MORE USEFUL DATA
               
                   ub::matrix<double>  rho_mat = ub::zero_matrix<double>(1,1);
            
                    
		    // evaluate AO Functions for all shells, NOW BLOCKWISE

                    // for each significant atom for this grid point
                    for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                                           
                     
                        // for each shell in this atom
                        for ( unsigned ishell = 0 ; ishell < _atomshells[rowatom].size() ; ishell++ ){
                            //boost::timer::cpu_times tstartshells = cpu_t.elapsed();
                            AOBasis::AOShellIterator _row = _atomshells[rowatom][ishell];
                            // for density, fill sub-part of AOatgrid
                           
                            ub::matrix_range< ub::matrix<double> > _AOgridsub = ub::subrange(AOgrid, 0, 1, (*_row)->getStartIndex(), (*_row)->getStartIndex()+(*_row)->getNumFunc());
                         
                            (*_row)->EvalAOspace(_AOgridsub, _grid[i][j].grid_pos);
                           

                        }  // shell in atom
                    }
                       
                   for ( unsigned sigrow = 0; sigrow < _significant_atoms[i][j].size() ; sigrow++){
                    
                        // this atom
                        int rowatom = _significant_atoms[i][j][sigrow];
                        ub::matrix<double> _temp     = ub::zero_matrix<double>(1,_blocksize[rowatom]);
                                          
                        ub::matrix_range< ub::matrix<double> > _AOgridrow     = ub::subrange(    AOgrid, 0,1, _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);

                        // for each atom
                        
                        for ( unsigned sigcol = 0; sigcol < _significant_atoms[i][j].size() ; sigcol++){
                            int colatom = _significant_atoms[i][j][sigcol];
                            
                            
                            // get the already calculated AO values
                           
                            ub::matrix_range< ub::matrix<double> >     _AOgridcol = ub::subrange(    AOgrid, 0, 1, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom]);
                        
                            ub::matrix_range<const ub::matrix<double> > DMAT_here = ub::subrange( _density_matrix, _startIdx[colatom], _startIdx[colatom]+_blocksize[colatom], _startIdx[rowatom], _startIdx[rowatom]+_blocksize[rowatom]);
                             
                            _temp     += ub::prod( _AOgridcol, DMAT_here);
                            
                            
                        } //col shells
                        


                        rho_mat  += ub::prod(_temp, ub::trans( _AOgridrow) );
                                               
                    } // row shells 


                    _grid[i][j].grid_density  =rho_mat(0,0);
                    Density_thread[i_thread] += _grid[i][j].grid_weight * _grid[i][j].grid_density;


                } // j: for each point in atom grid
                }// each thread
            } // i: for each atom grid

             for ( int i_thread = 0 ; i_thread < nthreads; i_thread++ ){
                
                result += Density_thread[i_thread]; 
                }
            density_set=true;
            return result;
         }
        

  
        double NumericalIntegration::StupidIntegrate(std::vector<double>& _data){
            
            
            double integral = 0.0;
            int _i_point = 0;
            for ( unsigned i = 0 ; i < _grid.size(); i++){
                for ( unsigned j = 0 ; j < _grid[i].size(); j++){

                    
                    integral += _data[_i_point] * _grid[i][j].grid_weight;
                    
                    _i_point++;

                }
            }
            
            return integral;
            
        }          
        
        std::vector<const vec *> NumericalIntegration::getGridpoints(){
            
            std::vector<const vec *> gridpoints;
            
            
            for ( unsigned i = 0 ; i < _grid.size(); i++){
                for ( unsigned j = 0 ; j < _grid[i].size(); j++){
                    gridpoints.push_back(&_grid[i][j].grid_pos);
                   
               }
                }
            return gridpoints;
        }
        
        
               
        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<ctp::QMAtom*> _atoms,AOBasis* basis) {
            _basis=basis;
            
            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers _grids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, _grids); // this checks out 1:1 with NWChem results! AWESOME

     
           map<string, GridContainers::radial_grid>::iterator it;

            LebedevGrid _sphericalgrid;
         
            for (it = _grids._radial_grids.begin(); it != _grids._radial_grids.end(); ++it) {
               _sphericalgrid.getSphericalGrid(_atoms, type, _grids);
       
            }

            
            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"
            
            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                double x_a = (*ait)->x * tools::conv::ang2bohr;
                double y_a = (*ait)->y * tools::conv::ang2bohr;
                double z_a = (*ait)->z * tools::conv::ang2bohr;
                int j = 0;
                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates in Bohr
                    double x_b = (*bit)->x * tools::conv::ang2bohr;
                    double y_b = (*bit)->y * tools::conv::ang2bohr;
                    double z_b = (*bit)->z * tools::conv::ang2bohr;

                    Rij.push_back(1.0 / sqrt((x_a - x_b)*(x_a - x_b) + (y_a - y_b)*(y_a - y_b) + (z_a - z_b)*(z_a - z_b)));


                    
                                        
                    j++;
                } // atoms
                Rij.push_back(0.0); // self-distance again
                i++;
            } // atoms
            


            int i_atom = 0;
            _totalgridsize = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                std::vector< GridContainers::integration_grid > _atomgrid;
                const vec atomA_pos =vec((*ait)->x * tools::conv::ang2bohr,(*ait)->y * tools::conv::ang2bohr,(*ait)->z * tools::conv::ang2bohr);
             
                string name = (*ait)->type;

                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = _grids._radial_grids.at(name);

                
                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = _grids._spherical_grids.at(name);

                // maximum order (= number of points) in spherical integration grid
                int maxorder = _sphericalgrid.Type2MaxOrder(name,type);
                int maxindex = _sphericalgrid.getIndexFromOrder(maxorder);

                // for pruning of integration grid, get interval boundaries for this element
                std::vector<double> PruningIntervals = _radialgrid.getPruningIntervals( name );
              //  cout << " Pruning Intervals: " << PruningIntervals[0] << " " << PruningIntervals[1] << " " << PruningIntervals[2] << " " << PruningIntervals[3] << endl;
                
                int current_order = 0;
                // get spherical grid
                std::vector<double> _theta;
                std::vector<double> _phi;
                std::vector<double> _weight;

                // for each radial value
                for (unsigned _i_rad = 0; _i_rad < _radial_grid.radius.size(); _i_rad++) {
                    double r = _radial_grid.radius[_i_rad];
                    int order;
                    // which Lebedev order for this point?
                    if ( maxindex == 1 ) {
                        // smallest possible grid anyway, nothing to do
                        order = maxorder;
                    } else if ( maxindex == 2 ) {
                        // only three intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(1);//1;
                        } else if ( ( r >= PruningIntervals[0] ) && ( r < PruningIntervals[3] )   ){
                            order = _sphericalgrid.getOrderFromIndex(2);
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(1);
                        } // maxorder == 2
                    } else {
                        // five intervals
                        if ( r < PruningIntervals[0] ) {
                            order = _sphericalgrid.getOrderFromIndex(int(2));
                        } else if ( ( r >= PruningIntervals[0]) && ( r < PruningIntervals[1] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(4);
                        } else if ( ( r >= PruningIntervals[1]) && ( r < PruningIntervals[2] ) ) {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1, 4));
                        } else if ( (r >= PruningIntervals[2]) && ( r < PruningIntervals[3] ) ) {
                            order = maxorder;
                        } else {
                            order = _sphericalgrid.getOrderFromIndex(max(maxindex-1,1));
                        }
                    }                        


                    
                    // get new spherical grid, if order changed
                    if ( order != current_order ){
                        _theta.clear();
                        _phi.clear();
                        _weight.clear();
                        
                        _sphericalgrid.getUnitSphereGrid(order,_theta,_phi,_weight);
                        current_order = order;
                    }
                    
                    // for each (theta,phi)
                    // for (int _i_sph = 0; _i_sph < _spherical_grid.phi.size(); _i_sph++) {

                    for (unsigned _i_sph = 0; _i_sph < _phi.size(); _i_sph++) {

                        double p   = _phi[_i_sph] * pi / 180.0; // back to rad
                        double t   = _theta[_i_sph] * pi / 180.0; // back to rad
                        double ws  = _weight[_i_sph];

                        const vec s = vec(sin(p) * cos(t), sin(p) * sin(t),cos(p));
                     


                        GridContainers::integration_grid _gridpoint;
                        _gridpoint.grid_pos = atomA_pos+r*s;

                        _gridpoint.grid_weight = _radial_grid.weight[_i_rad] * ws;

                        _atomgrid.push_back(_gridpoint);


                    } // spherical gridpoints
                } // radial gridpoint


                // get all distances from grid points to centers
                std::vector< std::vector<double> > rq;
                // for each center
                for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {
                    // get center coordinates
                   const vec atom_pos = vec((*bit)->x * tools::conv::ang2bohr,(*bit)->y * tools::conv::ang2bohr,(*bit)->z * tools::conv::ang2bohr);


                    std::vector<double> temp;
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {

                        temp.push_back(abs(git->grid_pos-atom_pos));

                    } // gridpoint of _atomgrid
                    rq.push_back(temp); // rq[center][gridpoint]

                } // centers
                // cout << " Calculated all gridpoint distances to centers for " << i_atom << endl;
                
                // find nearest-neighbor of this atom
                double distNN = 1e10;

                vector< ctp::QMAtom* > ::iterator NNit;
                //int i_NN;

                // now check all other centers
                int i_b =0;
                for (bit = _atoms.begin(); bit != _atoms.end(); ++bit) {

                    if (bit != ait) {
                        // get center coordinates
                       
                        const vec atomB_pos=vec((*bit)->x * tools::conv::ang2bohr,(*bit)->y * tools::conv::ang2bohr,(*bit)->z * tools::conv::ang2bohr);
                        double distSQ = (atomA_pos-atomB_pos)*(atomA_pos-atomB_pos);

                        // update NN distance and iterator
                        if ( distSQ < distNN ) {
                            distNN = distSQ;
                            NNit = bit;
                            //i_NN = i_b;
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers
              
                for ( unsigned i_grid = 0; i_grid < _atomgrid.size() ; i_grid++){
                    //cout << " modifying point " << i_grid << endl;
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( i_grid, _atoms.size(),rq);
                    //cout << " partition for gridpoint " << i_grid << endl;
                    // check weight sum
                    double wsum = 0.0;
                    for (unsigned i =0 ; i < _p.size(); i++ ){
                        wsum += _p[i];
                    }
                    //cout << " sum of partition weights " << wsum << endl;
                    if ( wsum != 0.0 ){
                        
                        // update the weight of this grid point
                        _atomgrid[i_grid].grid_weight = _atomgrid[i_grid].grid_weight * _p[i_atom]/wsum;
                        //cout << " adjusting gridpoint weight "  << endl;
                    } else {
                        
                       cerr << "\nSum of partition weights of grid point " << i_grid << " of atom " << i_atom << " is zero! ";
                       throw std::runtime_error("\nThis should never happen!"); 
                        
                    }
                    

                } // partition weight for each gridpoint

                // now remove points from the grid with negligible weights
                
                for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end();) {
                    if (git->grid_weight < 1e-13 ) {
                        git = _atomgrid.erase(git);
                    } else {
                        ++git;
                    }
                }
                
              
                
               // cout << " Total size of integration grid for atom: " << i_atom << " : " << _atomgrid.size() << " from " << fullsize << endl;

                _totalgridsize += _atomgrid.size() ;
                _grid.push_back(_atomgrid);
                i_atom++;
            } // atoms

            
           
            
            FindsignificantAtoms();
            return;
        }
    
        std::vector<double> NumericalIntegration::SSWpartition(int igrid, int ncenters, std::vector< std::vector<double> >& rq){
            const double ass = 0.725;
            // initialize partition vector to 1.0
            std::vector<double> p(ncenters,1.0);
            
            const double tol_scr = 1e-10;
            const double leps    = 1e-6; 
            // go through centers
            for ( int i = 1; i < ncenters; i++ ){
                
                int ij = i*(i+1)/2 -1; // indexing magic
                double rag = rq[i][igrid] ;
                
                // through all other centers (one-directional)
                for (int j = 0; j < i; j++ ){
                    
                    ij++;
                    if ( ( std::abs(p[i]) > tol_scr  ) || ( std::abs(p[j]) > tol_scr  ) ){
                        
                      
                        
                        double mu = ( rag - rq[j][igrid] )*Rij[ij]; 
                        if ( mu > ass ) {
                            p[i] = 0.0;
                        } else if ( mu < -ass ) {
                            p[j] = 0.0;
                        } else {
                            
                            double sk;
                            if (std::abs(mu) < leps ) {
                                sk = -1.88603178008*mu + 0.5;
                            } else {
                                sk = erf1c(mu); 
                            }
                            if ( mu > 0.0 ) sk = 1.0 - sk;
                            p[j] = p[j] * sk;
                            p[i] = p[i] * (1.0-sk);
                                                
                        }   
                    }  
                }

            }
            
            return p;
        }

        double NumericalIntegration::erf1c(double x){
             
            const static double alpha_erf1=1.0/0.30;
            return 0.5*erfcc((x/(1.0-x*x))*alpha_erf1);              
        }
              
        double NumericalIntegration::erfcc(double x){
            
            double tau = 1.0/(1.0+0.5*std::abs(x));
            
            return tau*exp(-x*x-1.26551223 + 1.00002368*tau + 0.37409196*tau*tau 
            + 0.09678418*pow(tau,3) - 0.18628806*pow(tau,4) + 0.27886807*pow(tau,5) 
            -1.13520398*pow(tau,6) + 1.48851587*pow(tau,7)  -0.82215223*pow(tau,8) 
            + 0.17087277*pow(tau,9));   
        }
                                                                                                
    }
}
