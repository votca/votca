/*
 *            Copyright 2016 The MUSCET Development Team
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


#include <votca/xtp/gyration.h>
#include <boost/format.hpp>
#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/orbitals.h>
//#include <votca/xtp/units.h>
#include <votca/tools/linalg.h>

namespace votca { namespace xtp {

void Density2Gyration::Initialize(Property* options) {
    string key = Identify();

    _state    = options->get(key + ".state").as<string> (); 
    _state_no = options->get(key + ".statenumber").as<int> ();
    _spin     = options->get(key + ".spin").as<string> ();
    if ( options->exists(key+".ecp")) {
       _use_ecp=options->get(key + ".ecp").as<bool> ();
    }

    _integrationmethod     = options->get(key + ".integrationmethod").as<string> ();
  
    if (!(_integrationmethod=="numeric" || _integrationmethod=="analytic")){
        std::runtime_error("Method not recognized. Only numeric and analytic available");
    }
    if ( options->exists(key+".gridsize")) {
         _gridsize = options->get(key+".gridsize").as<string>();
         }
    else _gridsize="medium";
    if ( options->exists(key+".openmp")) {
         _openmp_threads = options->get(key+".openmp").as<int>();
         }
    else _openmp_threads=0;
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

  
   
}



void Density2Gyration::AnalyzeDensity( Orbitals & _orbitals ){
    int threads=1;
#ifdef _OPENMP
            if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads); 
            threads=omp_get_max_threads();
#endif
   CTP_LOG(ctp::logDEBUG, *_log) << "===== Running on "<< threads << " threads ===== " << flush;

        vector< ctp::QMAtom* > Atomlist =_orbitals.QMAtoms();
        std::vector< ctp::QMAtom* >::iterator at;
        for (at=Atomlist.begin();at<Atomlist.end();++at){
            ctp::QMAtom * atom=new ctp::QMAtom(*(*at));
            _Atomlist.push_back(atom);
        }
        ub::matrix<double> DMAT_tot;
        BasisSet bs;
        bs.LoadBasisSet(_orbitals.getDFTbasis());
        AOBasis basis;
        basis.AOBasisFill(&bs, _Atomlist );
        
        // Analyze geometry
        AnalyzeGeometry( _Atomlist );
        
       
        
        std::vector<ub::matrix<double> > DMAT;

        //basis.ReorderMOs(_orbitals.MOCoefficients(), _orbitals.getQMpackage(), "votca" );  
        
        if(_state=="transition"){
                DMAT_tot=_orbitals.TransitionDensityMatrix(_spin, _state_no-1); 
        }
        else if (_state=="ground" || _state=="excited" || _state=="exciton" ){
             CTP_LOG(ctp::logDEBUG, *_log) << "Calculating density matrix:        " << _state << " No. " << _state_no << flush;
            
        
           
            ub::matrix<double> DMATGS=_orbitals.DensityMatrixGroundState();
            DMAT_tot=DMATGS;
            if ( _state_no > 0 && ( _state=="excited" || _state=="exciton" ) ){
               
                DMAT = _orbitals.DensityMatrixExcitedState( _spin, _state_no-1);
                
                if (_state == "excited" ){ 
                    DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1];
                }
               
            }            
	   // Ground state + hole_contribution + electron contribution
	}
        else throw std::runtime_error("State entry not recognized");

        
        if (_integrationmethod=="numeric")  {

            // setup numerical integration grid
            NumericalIntegration numway;
            numway.GridSetup(_gridsize,&bs,_Atomlist,&basis);
            
            
            if ( _state=="ground" || _state=="excited") {
                //LOG(logDEBUG, *_log) << TimeStamp() << " Calculate Densities at Numerical Grid with gridsize "<< _gridsize  << flush; 
                ub::vector<double> _analysis=numway.IntegrateGyrationTensor(DMAT_tot);

                // convert to eigenframe
                ub::vector<double> _gyration_tensor_diagonal;
                ub::matrix<double> _gyration_tensor_eigenframe;
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush; 
                Convert2Eigenframe( _analysis, _gyration_tensor_diagonal, _gyration_tensor_eigenframe );
            
                // determine quaternion for rotation of xyz to EF
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating Quaternion " << flush; 
                ub::vector<double> _quaternion = get_quaternion( _gyration_tensor_eigenframe );

                // report results
                CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush; 
                ReportAnalysis( _state, _analysis, _gyration_tensor_diagonal, _gyration_tensor_eigenframe  );
                
            
            } else if ( _state == "exciton" ){
                // hole density first
                ub::vector<double> _analysis_hole=numway.IntegrateGyrationTensor(DMAT[0]);

                // convert to eigenframe
                ub::vector<double> _gyration_tensor_diagonal_hole;
                ub::matrix<double> _gyration_tensor_eigenframe_hole;
                Convert2Eigenframe( _analysis_hole, _gyration_tensor_diagonal_hole, _gyration_tensor_eigenframe_hole );
            
                // determine quaternion for rotation of xyz to EF
                ub::vector<double> _quaternion_hole = get_quaternion( _gyration_tensor_eigenframe_hole );

                // report results
                ReportAnalysis( "hole", _analysis_hole, _gyration_tensor_diagonal_hole, _gyration_tensor_eigenframe_hole  );

                // electron density
                ub::vector<double> _analysis_electron=numway.IntegrateGyrationTensor(DMAT[1]);

                // convert to eigenframe
                ub::vector<double> _gyration_tensor_diagonal_electron;
                ub::matrix<double> _gyration_tensor_eigenframe_electron;
                Convert2Eigenframe( _analysis_electron, _gyration_tensor_diagonal_electron, _gyration_tensor_eigenframe_electron );
            
                // determine quaternion for rotation of xyz to EF
                 ub::vector<double> _quaternion_electron = get_quaternion( _gyration_tensor_eigenframe_electron );
                
                
                // report results
                ReportAnalysis( "electron", _analysis_electron, _gyration_tensor_diagonal_electron,  _gyration_tensor_eigenframe_electron  );
                
            }

        }
          else if (_integrationmethod=="analytic") {
              //esp.Fit2Density_analytic(_Atomlist,DMAT_tot,basis);
        }
        }


    void Density2Gyration::AnalyzeGeometry(vector<ctp::QMAtom*> _atoms){
    
        Elements _elements; 
        ub::vector<double> _analysis = ub::zero_vector<double>(10);
        std::vector< ctp::QMAtom* >::iterator at;
        for (at=_atoms.begin();at<_atoms.end();++at){
            
            double m = _elements.getMass( (*at)->type );
            double x = (*at)->x ;
            double y = (*at)->y ;
            double z = (*at)->z ;
            _analysis(0) += m   ;
            _analysis(1) += m*x ;
            _analysis(2) += m*y ;
            _analysis(3) += m*z ;
            
            _analysis(4) += m*x*x;
            _analysis(5) += m*x*y;
            _analysis(6) += m*x*z;
            _analysis(7) += m*y*y;
            _analysis(8) += m*y*z;
            _analysis(9) += m*z*z;
            
        }
        
        // normalize
        for ( unsigned i =1 ; i < 4; i++){
            _analysis(i) = _analysis(i)/_analysis(0)/tools::conv::bohr2ang;
        }
                // normalize
        for ( unsigned i =4 ; i < _analysis.size(); i++){
            _analysis(i) = _analysis(i)/_analysis(0)/tools::conv::bohr2ang/tools::conv::bohr2ang;
        }
        
        
        
        // gyration tensor
        _analysis(4) -= _analysis(1)*_analysis(1);
        _analysis(5) -= _analysis(1)*_analysis(2);
        _analysis(6) -= _analysis(1)*_analysis(3);
        _analysis(7) -= _analysis(2)*_analysis(2);
        _analysis(8) -= _analysis(2)*_analysis(3);
        _analysis(9) -= _analysis(3)*_analysis(3);

        
        // convert to eigenframe
        ub::vector<double> _gyration_tensor_diagonal;
        ub::matrix<double> _gyration_tensor_eigenframe;
        Convert2Eigenframe( _analysis, _gyration_tensor_diagonal, _gyration_tensor_eigenframe );

        cout << "\n";
        
        for ( int i =0 ; i< 3; i++){
            
            cout << _gyration_tensor_eigenframe(0,i) << " " << _gyration_tensor_eigenframe(1,i) << "  " <<  _gyration_tensor_eigenframe(2,i) << "\n" << endl;  
            
        }
        
        /*
        // test for right-handedness
        double check_x = _gyration_tensor_eigenframe(1,0)* _gyration_tensor_eigenframe(2,1) -  _gyration_tensor_eigenframe(2,0)* _gyration_tensor_eigenframe(1,1) ;
        double check_y = _gyration_tensor_eigenframe(2,0)* _gyration_tensor_eigenframe(0,1) -  _gyration_tensor_eigenframe(0,0)* _gyration_tensor_eigenframe(2,1) ;
        double check_z = _gyration_tensor_eigenframe(0,0)* _gyration_tensor_eigenframe(1,1) -  _gyration_tensor_eigenframe(1,0)* _gyration_tensor_eigenframe(0,1) ;
        
        cout << " check x" << check_x << " vs " << _gyration_tensor_eigenframe(0,2) << "\n" <<endl;
        cout << " check y" << check_y << " vs " << _gyration_tensor_eigenframe(1,2) << "\n" <<endl;
        cout << " check z" << check_z << " vs " << _gyration_tensor_eigenframe(2,2) << "\n" <<endl;
        */
        // determine quaternion for rotation of xyz to EF
        ub::vector<double> _quaternion = get_quaternion( _gyration_tensor_eigenframe );

        // report results
        ReportAnalysis( "geometry", _analysis, _gyration_tensor_diagonal, _gyration_tensor_eigenframe  );
        
        
        
    
    
    }



    void Density2Gyration::ReportAnalysis(string label, ub::vector<double> _tensor_elements, ub::vector<double> _tensor_diagonal, ub::matrix<double> _tensor_frame){

        
            CTP_LOG(ctp::logINFO, *_log) << "---------------- " << label << " ----------------" << flush;
            CTP_LOG(ctp::logINFO, *_log) << (boost::format("  Norm               = %1$9.4f ") % (_tensor_elements(0)) ) << flush;
            
            //LOG(logINFO,*_pLog) << (format("  Level = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = %4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") % (_i+_qpmin+1) % _dft_energies( _i + _qpmin ) % _vxc(_i,_i) % _sigma_x(_i,_i) % _sigma_c(_i,_i) % _qp_energies(_i + _qpmin ) ).str() << flush;
            
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid x         = %1$9.4f Ang") % (_tensor_elements(1)*tools::conv::bohr2ang) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid y         = %1$9.4f Ang") % (_tensor_elements(2)*tools::conv::bohr2ang) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid y         = %1$9.4f Ang") % (_tensor_elements(3)*tools::conv::bohr2ang) ) << flush;
            
            double RA2 = tools::conv::bohr2ang  *tools::conv::bohr2ang;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xx = %1$9.4f Ang^2") % (_tensor_elements(4)*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xy = %1$9.4f Ang^2") % (_tensor_elements(5)*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xz = %1$9.4f Ang^2") % (_tensor_elements(6)*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor yy = %1$9.4f Ang^2") % (_tensor_elements(7)*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor yz = %1$9.4f Ang^2") % (_tensor_elements(8)*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor zz = %1$9.4f Ang^2") % (_tensor_elements(9)*RA2) ) << flush;
            
            
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D1 = %1$9.4f Ang^2") % (_tensor_diagonal(0)*RA2) ) << flush;     
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D2 = %1$9.4f Ang^2") % (_tensor_diagonal(1)*RA2) ) << flush;     
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D3 = %1$9.4f Ang^2") % (_tensor_diagonal(2)*RA2) ) << flush;   

            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Radius of Gyration = %1$9.4f Ang") % (std::sqrt(_tensor_diagonal(0) + _tensor_diagonal(1) + _tensor_diagonal(2))*tools::conv::bohr2ang )) << flush;  
            
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 1 = %1$9.4f ") % (_tensor_frame(0,0)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 2 = %1$9.4f ") % (_tensor_frame(1,0)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 3 = %1$9.4f ") % (_tensor_frame(2,0)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 1 = %1$9.4f ") % (_tensor_frame(0,1)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 2 = %1$9.4f ") % (_tensor_frame(1,1)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 3 = %1$9.4f ") % (_tensor_frame(2,1)) ) << flush; 
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 1 = %1$9.4f ") % (_tensor_frame(0,2)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 2 = %1$9.4f ") % (_tensor_frame(1,2)) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 3 = %1$9.4f ") % (_tensor_frame(2,2)) ) << flush;             

    
    }

    ub::vector<double> Density2Gyration::get_quaternion(ub::matrix<double>& eigenframe){
        
        // second coordinate system is assumed to be Cartesian
     
        ub::matrix<double> M=ub::zero_matrix<double>(3,3);

        // add all the outer products, if second coordinate system was not Cartesian
        // is anyways, so M=eigenframe
        //M=ub::trans(eigenframe);
        M=eigenframe;
       
        ub::matrix<double> N = ub::zero_matrix<double>(4,4);
          
        N(0,0) =  M(0,0) + M(1,1) + M(2,2); // N11=float(M[0][:,0]+M[1][:,1]+M[2][:,2]);
        N(1,1) =  M(0,0) - M(1,1) - M(2,2); // N22=float(M[0][:,0]-M[1][:,1]-M[2][:,2])
        N(2,2) = -M(0,0) + M(1,1) - M(2,2); // N33=float(-M[0][:,0]+M[1][:,1]-M[2][:,2])
        N(3,3) = -M(0,0) - M(1,1) + M(2,2); // N44=float(-M[0][:,0]-M[1][:,1]+M[2][:,2])
        N(0,1) =  M(1,2) - M(2,1);          // N12=float(M[1][:,2]-M[2][:,1])
        N(0,2) =  M(2,0) - M(0,2);          // N13=float(M[2][:,0]-M[0][:,2])
        N(0,3) =  M(0,1) - M(1,0);          // N14=float(M[0][:,1]-M[1][:,0])
        N(1,0) =  N(0,1);                   // N21=float(N12)
        N(1,2) =  M(0,1) + M(1,0);          // N23=float(M[0][:,1]+M[1][:,0])
        N(1,3) =  M(2,0) + M(0,2);          // N24=float(M[2][:,0]+M[0][:,2])
        N(2,0) =  N(0,2);                   // N31=float(N13)
        N(2,1) =  N(1,2);                   // N32=float(N23)
        N(2,3) =  M(1,2) + M(2,1);          // N34=float(M[1][:,2]+M[2][:,1])
        N(3,0) =  N(0,3);                   // N41=float(N14)
        N(3,1) =  N(1,3);                   // N42=float(N24)
        N(3,2) =  N(2,3);                   // N43=float(N34)

        ub::vector<double> eigenvalues;
        ub::matrix<double> eigenvectors;
        linalg_eigenvalues( N, eigenvalues, eigenvectors );
        
        // find max eigenvalue
        int index=0;
        double maxev = eigenvalues(0);

        for (unsigned i = 1; i < eigenvalues.size(); i++ ) {
            if ( eigenvalues(i) > maxev ) {
                maxev = eigenvalues(i);
                index = i;
            }
        }
        
        // return this vector
        ub::vector<double> quaternion = ub::zero_vector<double>(4);
        quaternion(0)= eigenvectors(0,index);
        quaternion(1)= eigenvectors(1,index);
        quaternion(2)= eigenvectors(2,index);
        quaternion(3)= eigenvectors(3,index);
        
        return quaternion;
        
    }
    
    
    
    
        void Density2Gyration::Convert2Eigenframe(ub::vector<double> V, ub::vector<double> &_gyration_tensor_diagonal, ub::matrix<double> &_gyration_tensor_eigenframe   ){
            
                 ub::matrix<double> _gyration_tensor = ub::zero_matrix<double>(3,3);
                _gyration_tensor(0,0) = V(4);
                _gyration_tensor(1,0) = V(5);
                _gyration_tensor(2,0) = V(6);
                _gyration_tensor(0,1) = V(5);
                _gyration_tensor(1,1) = V(7);
                _gyration_tensor(2,1) = V(8);
                _gyration_tensor(0,2) = V(6);
                _gyration_tensor(1,2) = V(8);
                _gyration_tensor(2,2) = V(9);
            
                linalg_eigenvalues( _gyration_tensor, _gyration_tensor_diagonal, _gyration_tensor_eigenframe );
            
        }




}}


