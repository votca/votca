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
#include <numeric>

#include <votca/xtp/aomatrix.h>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>

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
                throw std::runtime_error("Specify at least one functional");
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
        
        ub::matrix<double> NumericalIntegration::IntegrateExternalPotential(const std::vector<double>& Potentialvalues){
            
            ub::matrix<double> ExternalMat = ub::zero_matrix<double>(_basis->AOBasisSize(), _basis->AOBasisSize());
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               std::vector<ub::matrix<double> >vex_thread;
               std::vector<double> Exc_thread=std::vector<double>(nthreads,0.0);
               for(unsigned i=0;i<nthreads;++i){
                   ub::matrix<double> Vex_thread=ub::zero_matrix<double>(ExternalMat.size1());
                   vex_thread.push_back(Vex_thread);
               }
               
               
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
            for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {
                
                
                GridBox& box = _grid_boxes[i];
                
               
                
                ub::matrix<double> Vex_here=ub::zero_matrix<double>(box.Matrixsize());
                const std::vector<tools::vec>& points=box.getGridPoints();
                const std::vector<double>& weights=box.getGridWeights();
                
                ub::range one=ub::range(0,1);
                
                ub::matrix<double> _temp     = ub::zero_matrix<double>(1,box.Matrixsize());
                
                ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                
                
                
                //iterate over gridpoints
                for(unsigned p=0;p<box.size();p++){
                    ao=ub::zero_matrix<double>(1,box.Matrixsize());
                    const std::vector<ub::range>& aoranges=box.getAOranges();
                    const std::vector<const AOShell* > shells=box.getShells();
                    for(unsigned j=0;j<box.Shellsize();++j){
                        const AOShell* shell=shells[j];
                        ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao,one,aoranges[j]);
                        shell->EvalAOspace(aoshell,points[p]);
                    }

                    double weight=weights[p];
                    ub::matrix<double> _addEX = weight*Potentialvalues[box.getIndexoffirstgridpoint()+p]*ao ;
                    
                    Vex_here+=ub::prod( ub::trans(_addEX), ao);
                }
                
                
                box.AddtoBigMatrix(vex_thread[thread],Vex_here);
                
            }
            }   
            for(unsigned i=0;i<nthreads;++i){
                ExternalMat+=vex_thread[i];
                
               }   
         
         
            
            ExternalMat += ub::trans(ExternalMat);
            return ExternalMat;

        }
        
        
        void NumericalIntegration::setXCfunctional(const string _functional){
            
            Vxc_Functionals map;
            std::vector<string> strs;           
            boost::split(strs, _functional, boost::is_any_of(" "));
            xfunc_id = 0;
            
#ifdef LIBXC
            _use_votca = false;
            _use_separate = false;
            cfunc_id = 0;

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
            if(_use_votca){
                cout<<"Warning: VOTCA_PBE does give correct Vxc but incorrect E_xc"<<endl;
            }
            setXC=true;
            return;
        }

        
        
        
        void NumericalIntegration::EvaluateXC(const double rho,const ub::matrix<double>& grad_rho,double& f_xc, double& df_drho, double& df_dsigma){
            
                              
 #ifdef LIBXC                   
                    if (_use_votca) {
#endif                  
                        _xc.getXC(xfunc_id, rho, grad_rho(0,0), grad_rho(0,1), grad_rho(0,2), f_xc, df_drho, df_dsigma);
#ifdef LIBXC
                    }                        // evaluate via LIBXC, if compiled, otherwise, go via own implementation

                    else {
                     

                        double sigma = ub::prod(grad_rho,ub::trans(grad_rho))(0,0);

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
            
            return;
        }
            
        double NumericalIntegration::IntegratePotential(const vec& rvector){
            
            double result = 0.0;
            assert(density_set && "Density not calculated");
          
            for (unsigned i = 0; i < _grid_boxes.size(); i++) {

                const std::vector<tools::vec>& points = _grid_boxes[i].getGridPoints();
                const std::vector<double>& weights = _grid_boxes[i].getGridWeights();
                const std::vector<double>& densities = _grid_boxes[i].getGridDensities();
                for (unsigned j = 0; j < points.size(); j++) {
                    double dist = abs(points[j] - rvector);
                    result -= weights[j] * densities[j] / dist;
                }
            }

            return result;   
        }
               
        
        
        
        
        void NumericalIntegration::SortGridpointsintoBlocks(std::vector< std::vector< GridContainers::integration_grid > >& grid){
            const double boxsize=2.5;
            
            std::vector< std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > >  boxes;
            
            tools::vec min=vec(std::numeric_limits<double>::max());
            tools::vec max=vec(std::numeric_limits<double>::min());
                   
            for ( unsigned i = 0 ; i < grid.size(); i++){
                for ( unsigned j = 0 ; j < grid[i].size(); j++){
                    const tools::vec& pos= grid[i][j].grid_pos;
                    if(pos.getX()>max.getX()){
                        max.x()=pos.getX();
                    }
                    else if(pos.getX()<min.getX()){
                        min.x()=pos.getX();
                    }
                    if(pos.getY()>max.getY()){
                        max.y()=pos.getY();
                    }
                    else if(pos.getY()<min.getY()){
                        min.y()=pos.getY();
                    }
                    if(pos.getZ()>max.getZ()){
                        max.z()=pos.getZ();
                    }
                    else if(pos.getZ()<min.getZ()){
                        min.z()=pos.getZ();
                        }
                    }
                }
            
            vec molextension=(max-min);
            vec numberofboxes=molextension/boxsize;
            vec roundednumofbox=vec(std::ceil(numberofboxes.getX()),std::ceil(numberofboxes.getY()),std::ceil(numberofboxes.getZ()));

            
            //creating temparray
            for (unsigned i=0;i<unsigned(roundednumofbox.getX());i++){
                std::vector< std::vector< std::vector< GridContainers::integration_grid* > > > boxes_yz;
                for (unsigned j=0;j<unsigned(roundednumofbox.getY());j++){
                    std::vector< std::vector< GridContainers::integration_grid* > >  boxes_z;
                    for (unsigned k=0;k<unsigned(roundednumofbox.getZ());k++){
                        std::vector< GridContainers::integration_grid* >  box;
                        box.reserve(100);
                        boxes_z.push_back(box);
                    }
                    boxes_yz.push_back(boxes_z);
            }
                boxes.push_back(boxes_yz);
            }
            
             for ( auto & atomgrid : grid){
                for ( auto & gridpoint : atomgrid){
                    tools::vec pos= gridpoint.grid_pos-min;
                    tools::vec index=pos/boxsize;
                    int i_x=int(index.getX());
                    int i_y=int(index.getY());
                    int i_z=int(index.getZ());
                    boxes[i_x][i_y][i_z].push_back(&gridpoint);
                }
             }
            
            for ( auto& boxes_xy : boxes){
                for( auto& boxes_z : boxes_xy){
                    for ( auto& box : boxes_z){      
                        if( box.size()<1){
                            continue;
                        }
                        GridBox gridbox;
                        
                        for(const auto&point:box){
                            gridbox.addGridPoint(*point);
                        }
                        _grid_boxes.push_back(gridbox);
                    }
                }
            }
            
            return;
        }
        
        
        void NumericalIntegration::FindSignificantShells(){

            for (unsigned i=0;i<_grid_boxes.size();++i){
                GridBox & box=_grid_boxes[i];
                for (AOBasis::AOShellIterator _row = _basis->firstShell(); _row != _basis->lastShell(); _row++) {
                      AOShell* _store=(*_row);
                      const double decay=(*_row)->getMinDecay();
                      const tools::vec& shellpos=(*_row)->getPos();
                      
                      for(const auto& point : box.getGridPoints()){
                          tools::vec dist=shellpos-point;
                          double distsq=dist*dist;
                          // if contribution is smaller than -ln(1e-10), add atom to list
                        if ( (decay * distsq) < 20.7 ){
                            box.addShell(_store);
                            break;
                        }
                      }
                }
                //cout<<box.significant_shells.size()<<" "<<box.grid_pos.size()<<endl;
            }
            
             std::vector< GridBox > _grid_boxes_copy;
            
            int combined=0;
            std::vector<bool> Compared=std::vector<bool>(_grid_boxes.size(),false);
            for (unsigned i=0;i<_grid_boxes.size();i++){
                if(Compared[i]){continue;}
                GridBox box=_grid_boxes[i];
                if(box.Shellsize()<1){continue;}
                Compared[i]=true;
                for (unsigned j=i+1;j<_grid_boxes.size();j++){                   
                    if(GridBox::compareGridboxes(_grid_boxes[i],_grid_boxes[j])){
                        Compared[j]=true;
                        box.addGridBox(_grid_boxes[j]);
                        combined++;
                    }
                    
                }
                _grid_boxes_copy.push_back(box);
            }

         
            
            
            std::vector<unsigned> sizes;
            sizes.reserve(_grid_boxes_copy.size());
            for(auto& box: _grid_boxes_copy){
                sizes.push_back(box.size()*box.Matrixsize());
            }
           
            
            std::vector<unsigned> indexes=std::vector<unsigned>(sizes.size());
            std::iota(indexes.begin(), indexes.end(), 0);
            std::sort(indexes.begin(), indexes.end(),[&sizes](unsigned i1, unsigned i2) {return sizes[i1] > sizes[i2];});
            
             unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
            
            std::vector<unsigned> scores=std::vector<unsigned>(nthreads,0);
            std::vector< std::vector<unsigned> > indices;
            for (unsigned i=0;i<nthreads;++i){
                std::vector<unsigned> thread_box_indices;
                indices.push_back(thread_box_indices);
            }
        
            
            for(const auto index:indexes){
                unsigned thread=0;
                unsigned minimum= std::numeric_limits<unsigned>::max();
                for(unsigned i=0;i<scores.size();++i){
                    if(scores[i]<minimum){
                        minimum=scores[i];
                        thread=i;
                    }
                }
                indices[thread].push_back(index);
                scores[thread]+=sizes[index];   
            }           
            
            thread_start=std::vector<unsigned>(0);
            thread_stop=std::vector<unsigned>(0);
            unsigned start=0;
            unsigned stop=0;
            unsigned indexoffirstgridpoint=0;
             _grid_boxes.resize(0);
            for (const std::vector<unsigned>& thread_index:indices){
                thread_start.push_back(start);
                stop=start+thread_index.size();
                thread_stop.push_back(stop);
                start=stop;
                for(const unsigned index:thread_index){        
                        GridBox newbox=_grid_boxes_copy[index];
                        newbox.setIndexoffirstgridpoint(indexoffirstgridpoint);
                        indexoffirstgridpoint+=newbox.size();
                        newbox.PrepareForIntegration();
                        _grid_boxes.push_back(newbox);                 
                }
            }   
            return;
        }
        
        
        
        
        ub::matrix<double> NumericalIntegration::IntegrateVXC(const ub::matrix<double>& _density_matrix){
            ub::matrix<double> Vxc=ub::zero_matrix<double>(_density_matrix.size1());
            EXC = 0;
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               std::vector<ub::matrix<double> >vxc_thread;
               std::vector<double> Exc_thread=std::vector<double>(nthreads,0.0);
               for(unsigned i=0;i<nthreads;++i){
                   ub::matrix<double> Vxc_thread=ub::zero_matrix<double>(_density_matrix.size1());
                   vxc_thread.push_back(Vxc_thread);
               }
               
               
            #pragma omp parallel for
             for (unsigned thread=0;thread<nthreads;++thread){
                for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {
                
                double EXC_box=0.0;
                GridBox& box = _grid_boxes[i];
               
                const ub::matrix<double>  DMAT_here=box.ReadFromBigMatrix(_density_matrix);
                
                ub::matrix<double> Vxc_here=ub::zero_matrix<double>(DMAT_here.size1());
                const std::vector<tools::vec>& points=box.getGridPoints();
                const std::vector<double>& weights=box.getGridWeights();
                
                ub::range one=ub::range(0,1);
                ub::range three=ub::range(0,3);
                ub::matrix<double> _temp     = ub::zero_matrix<double>(1,box.Matrixsize());
                ub::matrix<double> _tempgrad = ub::zero_matrix<double>(3,box.Matrixsize());
                ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                ub::matrix<double> ao_grad=ub::matrix<double>(3,box.Matrixsize());
               
                //iterate over gridpoints
                for(unsigned p=0;p<box.size();p++){
                    ao=ub::zero_matrix<double>(1,box.Matrixsize());
                    ao_grad=ub::zero_matrix<double>(3,box.Matrixsize());
                    const std::vector<ub::range>& aoranges=box.getAOranges();
                    const std::vector<const AOShell* >& shells=box.getShells();
                   
                    for(unsigned j=0;j<box.Shellsize();++j){
                        const AOShell* shell=shells[j];
                       
                        ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao,one,aoranges[j]);
                       
                        ub::matrix_range< ub::matrix<double> > ao_grad_shell=ub::project(ao_grad,three,aoranges[j]);
                        
                        shell->EvalAOspace(aoshell,ao_grad_shell,points[p]);
                       
                    }
                    
                    _temp=ub::prod( ao, DMAT_here);
                   
                    _tempgrad=ub::prod(ao_grad,DMAT_here);
                    
                    double rho=ub::prod(_temp, ub::trans( ao) )(0,0);
                    
                    ub::matrix<double> rho_grad=ub::prod(_temp, ub::trans(ao_grad))+ub::prod(ao,ub::trans(_tempgrad));
                   
		    if ( rho < 1.e-15 ) continue; // skip the rest, if density is very small
                    
                    double f_xc;      // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r }
                    double df_drho;   // v_xc_rho(r) = df/drho
                    double df_dsigma; // df/dsigma ( df/dgrad(rho) = df/dsigma * dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
                    EvaluateXC( rho,rho_grad,f_xc, df_drho, df_dsigma);
                    
                    double weight=weights[p];
                    ub::matrix<double> _addXC = weight * df_drho * ao *0.5;
                    
                    _addXC+=  2.0*df_dsigma * weight * ub::prod(rho_grad,ao_grad);

                    // Exchange correlation energy
                    EXC_box += weight  * rho * f_xc;
                  
                    Vxc_here+=ub::prod( ub::trans(_addXC), ao);
                  
                }
                
                box.AddtoBigMatrix(vxc_thread[thread],Vxc_here);
              
                Exc_thread[thread]+=EXC_box;
                
            }
                
            }   
            for(unsigned i=0;i<nthreads;++i){
                Vxc+=vxc_thread[i];
                EXC+=Exc_thread[i];
               }   
            Vxc+=ub::trans(Vxc);
            
            return Vxc;
        }
        
        double NumericalIntegration::IntegrateDensity(const ub::matrix<double>& _density_matrix){
            
            double N = 0;
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               
               std::vector<double> N_thread=std::vector<double>(nthreads,0.0);
               
               
               
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
             for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {
                
                double N_box=0.0;
                GridBox& box = _grid_boxes[i];
                
                
                const ub::matrix<double>  DMAT_here=box.ReadFromBigMatrix(_density_matrix);
                
                ub::matrix<double> Vxc_here=ub::zero_matrix<double>(DMAT_here.size1());
                const std::vector<tools::vec>& points=box.getGridPoints();
                const std::vector<double>& weights=box.getGridWeights();
                
                ub::range one=ub::range(0,1);
                
                ub::matrix<double> _temp     = ub::zero_matrix<double>(1,box.Matrixsize());
               
                ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                
                box.prepareDensity();
                
                //iterate over gridpoints
                for(unsigned p=0;p<box.size();p++){
                    ao=ub::zero_matrix<double>(1,box.Matrixsize());
                   
                    const std::vector<ub::range>& aoranges=box.getAOranges();
                    const std::vector<const AOShell* > shells=box.getShells();
                    for(unsigned j=0;j<box.Shellsize();++j){
                        const AOShell* shell=shells[j];
                        ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao,one,aoranges[j]);
                        
                        shell->EvalAOspace(aoshell,points[p]);
                    }
                    
                    _temp=ub::prod( ao, DMAT_here);
                   
                    
                    
                    double rho=ub::prod(_temp, ub::trans( ao) )(0,0);
                    box.addDensity(rho);
                    N_box+=rho*weights[p];
                    
                }

                N_thread[thread]+=N_box;
                
            }
            }   
            for(unsigned i=0;i<nthreads;++i){
                N+=N_thread[i];
               }   
            density_set=true;
            return N;
        }
        
        
 ub::vector<double> NumericalIntegration::IntegrateGyrationTensor(const ub::matrix<double>& _density_matrix){
            
            double N = 0;
            double centroid_x = 0.0;
            double centroid_y = 0.0;
            double centroid_z = 0.0;
            double gyration_xx = 0.0;
            double gyration_xy = 0.0;
            double gyration_xz = 0.0;
            double gyration_yy = 0.0;
            double gyration_yz = 0.0;
            double gyration_zz = 0.0;
            ub::vector<double> result=ub::zero_vector<double>(10);
            
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               
               std::vector<double> N_thread=std::vector<double>(nthreads,0.0);

               // centroid
	       std::vector<double> centroid_x_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> centroid_y_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> centroid_z_thread=std::vector<double>(nthreads,0.0);

	       // gyration tensor
	       std::vector<double> gyration_xx_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> gyration_xy_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> gyration_xz_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> gyration_yy_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> gyration_yz_thread=std::vector<double>(nthreads,0.0);
	       std::vector<double> gyration_zz_thread=std::vector<double>(nthreads,0.0);
               
               
            #pragma omp parallel for
            for (unsigned thread=0;thread<nthreads;++thread){
             for (unsigned i = thread_start[thread]; i < thread_stop[thread]; ++i) {
                
                double N_box=0.0;
                double centroid_x_box=0.0;
                double centroid_y_box=0.0;
                double centroid_z_box=0.0;
                double gyration_xx_box=0.0;
                double gyration_xy_box=0.0;
                double gyration_xz_box=0.0;
                double gyration_yy_box=0.0;
                double gyration_yz_box=0.0;
                double gyration_zz_box=0.0;
                
                GridBox& box = _grid_boxes[i];
                
                
                const ub::matrix<double>  DMAT_here=box.ReadFromBigMatrix(_density_matrix);
                
                const std::vector<tools::vec>& points=box.getGridPoints();
                const std::vector<double>& weights=box.getGridWeights();
                
                ub::range one=ub::range(0,1);
                
                ub::matrix<double> _temp     = ub::zero_matrix<double>(1,box.Matrixsize());
               
                ub::matrix<double> ao=ub::matrix<double>(1,box.Matrixsize());
                
                box.prepareDensity();
                
                //iterate over gridpoints
                for(unsigned p=0;p<box.size();p++){
                    ao=ub::zero_matrix<double>(1,box.Matrixsize());
                   
                    const std::vector<ub::range>& aoranges=box.getAOranges();
                    const std::vector<const AOShell* > shells=box.getShells();
                    for(unsigned j=0;j<box.Shellsize();++j){
                        const AOShell* shell=shells[j];
                        ub::matrix_range< ub::matrix<double> > aoshell=ub::project(ao,one,aoranges[j]);
                        
                        shell->EvalAOspace(aoshell,points[p]);
                    }
                    
                    _temp=ub::prod( ao, DMAT_here);
                   
                    
                    
                    double rho=ub::prod(_temp, ub::trans( ao) )(0,0);
                    box.addDensity(rho);
                    N_box+=rho*weights[p];
                    centroid_x_box+=rho*weights[p]* points[p].getX();
                    centroid_y_box+=rho*weights[p]* points[p].getY();
                    centroid_z_box+=rho*weights[p]* points[p].getZ();
                    gyration_xx_box+=rho*weights[p]* points[p].getX()*points[p].getX();
                    gyration_xy_box+=rho*weights[p]* points[p].getX()*points[p].getY();
                    gyration_xz_box+=rho*weights[p]* points[p].getX()*points[p].getZ();
                    gyration_yy_box+=rho*weights[p]* points[p].getY()*points[p].getY();
                    gyration_yz_box+=rho*weights[p]* points[p].getY()*points[p].getZ();
                    gyration_zz_box+=rho*weights[p]* points[p].getZ()*points[p].getZ();
                    
                }

                N_thread[thread]+=N_box;
                centroid_x_thread[thread] += centroid_x_box;
                centroid_y_thread[thread] += centroid_y_box;
                centroid_z_thread[thread] += centroid_z_box;
                gyration_xx_thread[thread] += gyration_xx_box;
                gyration_xy_thread[thread] += gyration_xy_box;
                gyration_xz_thread[thread] += gyration_xz_box;
                gyration_yy_thread[thread] += gyration_yy_box;
                gyration_yz_thread[thread] += gyration_yz_box;
                gyration_zz_thread[thread] += gyration_zz_box;
                
            }
            }   
            for(unsigned i=0;i<nthreads;++i){
                N+=N_thread[i];
                centroid_x += centroid_x_thread[i];
                centroid_y += centroid_y_thread[i];
                centroid_z += centroid_z_thread[i];
                gyration_xx += gyration_xx_thread[i];
                gyration_xy += gyration_xy_thread[i];
                gyration_xz += gyration_xz_thread[i];
                gyration_yy += gyration_yy_thread[i];
                gyration_yz += gyration_yz_thread[i];
                gyration_zz += gyration_zz_thread[i];
               }   
            density_set=true;

            // Normalize
	    centroid_x = centroid_x/N;
	    centroid_y = centroid_y/N;
	    centroid_z = centroid_z/N;

            gyration_xx = gyration_xx/N;
	    gyration_xy = gyration_xy/N;
	    gyration_xz = gyration_xz/N;
	    gyration_yy = gyration_yy/N;
	    gyration_yz = gyration_yz/N;
	    gyration_zz = gyration_zz/N;

	      // Copy all elements to result vector
	      result(0) = N;
	      result(1) = centroid_x;
	      result(2) = centroid_y;
	      result(3) = centroid_z;
	      result(4) = gyration_xx - centroid_x*centroid_x;
	      result(5) = gyration_xy - centroid_x*centroid_y;
	      result(6) = gyration_xz - centroid_x*centroid_z;
	      result(7) = gyration_yy - centroid_y*centroid_y;
	      result(8) = gyration_yz - centroid_y*centroid_z;
	      result(9) = gyration_zz - centroid_z*centroid_z;
            

            return result;
        }
        
  
        
        
        std::vector<const vec *> NumericalIntegration::getGridpoints(){
            
            std::vector<const vec *> gridpoints;
            
            
            for ( unsigned i = 0 ; i < _grid_boxes.size(); i++){
                const std::vector<tools::vec>& points=_grid_boxes[i].getGridPoints();
                for ( unsigned j = 0 ; j < points.size(); j++){
                    gridpoints.push_back(&points[j]);
                   
               }
                }
            return gridpoints;
        }
        
        
               
        void NumericalIntegration::GridSetup(string type, BasisSet* bs, vector<ctp::QMAtom*> _atoms,AOBasis* basis) {
            _basis=basis;
            std::vector< std::vector< GridContainers::integration_grid > > grid;
            const double pi = boost::math::constants::pi<double>();
            // get GridContainer
            GridContainers initialgrids;

            // get radial grid per element
            EulerMaclaurinGrid _radialgrid;
            _radialgrid.getRadialGrid(bs, _atoms, type, initialgrids); // this checks out 1:1 with NWChem results! AWESOME

     
           map<string, GridContainers::radial_grid>::iterator it;

            LebedevGrid _sphericalgrid;
         
            for (it = initialgrids._radial_grids.begin(); it != initialgrids._radial_grids.end(); ++it) {
               _sphericalgrid.getSphericalGrid(_atoms, type, initialgrids);
       
            }

            
            // for the partitioning, we need all inter-center distances later, stored in one-directional list
            int ij = 0;
            Rij.push_back(0.0); // 1st center "self-distance"
            
            vector< ctp::QMAtom* > ::iterator ait;
            vector< ctp::QMAtom* > ::iterator bit;
            int i = 1;
            for (ait = _atoms.begin() + 1; ait != _atoms.end(); ++ait) {
                // get center coordinates in Bohr
                vec pos_a = (*ait)->getPos() * tools::conv::ang2bohr;
                
                int j = 0;
                for (bit = _atoms.begin(); bit != ait; ++bit) {
                    ij++;
                    // get center coordinates in Bohr
                    vec pos_b = (*bit)->getPos() * tools::conv::ang2bohr;
                   
                    Rij.push_back(1.0 / abs(pos_a-pos_b));
                                        
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
                const vec atomA_pos =(*ait)->getPos() * tools::conv::ang2bohr;
             
                string name = (*ait)->type;
                
                // get radial grid information for this atom type
                GridContainers::radial_grid _radial_grid = initialgrids._radial_grids.at(name);

                
                // get spherical grid information for this atom type
                GridContainers::spherical_grid _spherical_grid = initialgrids._spherical_grids.at(name);

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
                   const vec atom_pos = (*bit)->getPos() * tools::conv::ang2bohr;


                    std::vector<double> temp;
                    // for each gridpoint
                    for (std::vector<GridContainers::integration_grid >::iterator git = _atomgrid.begin(); git != _atomgrid.end(); ++git) {

                        temp.push_back(abs(git->grid_pos-atom_pos));

                    } // gridpoint of _atomgrid
                    rq.push_back(temp); 

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
                       
                        const vec atomB_pos=(*bit)->getPos() * tools::conv::ang2bohr;
                        double distSQ = (atomA_pos-atomB_pos)*(atomA_pos-atomB_pos);

                        // update NN distance and iterator
                        if ( distSQ < distNN ) {
                            distNN = distSQ;
                            NNit = bit;
                           
                        }

                    } // if ( ait != bit) 
                    i_b++;
                }// bit centers
                
                for ( unsigned i_grid = 0; i_grid < _atomgrid.size() ; i_grid++){
                    // call some shit called grid_ssw0 in NWChem
                    std::vector<double> _p = SSWpartition( i_grid, _atoms.size(),rq);
                 
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
                
                _totalgridsize += _atomgrid.size() ;

                grid.push_back(_atomgrid);
                
                i_atom++;
                
            } // atoms
            
            SortGridpointsintoBlocks(grid);
            FindSignificantShells();
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
