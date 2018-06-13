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
#include <votca/tools/elements.h>


using namespace votca::tools;

namespace votca { namespace xtp {

void Density2Gyration::Initialize( tools::Property* options) {
    string key = Identify();

    _state    = options->get(key + ".state").as<string> (); 
    _state_no = options->get(key + ".statenumber").as<int> ();
    _spin     = options->get(key + ".spin").as<string> ();
 
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

    void Density2Gyration::AnalyzeDensity(Orbitals & _orbitals) {
      int threads = 1;
#ifdef _OPENMP
      if (_openmp_threads > 0) omp_set_num_threads(_openmp_threads);
      threads = omp_get_max_threads();
#endif
      CTP_LOG(ctp::logDEBUG, *_log) << "===== Running on " << threads << " threads ===== " << flush;

      std::vector< QMAtom* > _Atomlist = _orbitals.QMAtoms();
      Eigen::MatrixXd DMAT_tot;
      BasisSet bs;
      bs.LoadBasisSet(_orbitals.getDFTbasis());
      AOBasis basis;
      basis.AOBasisFill(&bs, _Atomlist);
      // Analyze geometry
      AnalyzeGeometry(_Atomlist);
      std::vector<Eigen::MatrixXd > DMAT;

      if (_state == "transition") {
        DMAT_tot = _orbitals.TransitionDensityMatrix(_spin, _state_no - 1);
      } else if (_state == "ground" || _state == "excited" || _state == "exciton") {
        CTP_LOG(ctp::logDEBUG, *_log) << "Calculating density matrix:        " << _state << " No. " << _state_no << flush;
        Eigen::MatrixXd DMATGS = _orbitals.DensityMatrixGroundState();
        DMAT_tot = DMATGS;
        if (_state_no > 0 && (_state == "excited" || _state == "exciton")) {
          DMAT = _orbitals.DensityMatrixExcitedState(_spin, _state_no - 1);
          if (_state == "excited") {
            DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1];
          }
        }
        // Ground state + hole_contribution + electron contribution
      } else throw std::runtime_error("State entry not recognized");

      // setup numerical integration grid
      NumericalIntegration numway;
      numway.GridSetup(_gridsize, _Atomlist, &basis);

      if (_state == "ground" || _state == "excited") {
        //LOG(logDEBUG, *_log) << TimeStamp() << " Calculate Densities at Numerical Grid with gridsize "<< _gridsize  << flush; 
        Gyrationtensor gyro = numway.IntegrateGyrationTensor(DMAT_tot);
        tools::matrix::eigensystem_t system;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro.gyration.SolveEigensystem(system); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating Quaternion " << flush;
        //Eigen::Quaterniond _quaternion = get_quaternion( system );
        // report results
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis(_state, gyro,system);


      } else if (_state == "exciton") {
        // hole density first
        Gyrationtensor gyro_hole = numway.IntegrateGyrationTensor(DMAT[0]);
        tools::matrix::eigensystem_t system_h;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro_hole.gyration.SolveEigensystem(system_h); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating Quaternion " << flush;
        //Eigen::Quaterniond _quaternion_h = get_quaternion(system_h );
        // report results
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis("hole", gyro_hole,system_h);
        
        // electron density
        Gyrationtensor gyro_electron = numway.IntegrateGyrationTensor(DMAT[1]);
        tools::matrix::eigensystem_t system_e;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro_electron.gyration.SolveEigensystem(system_e); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating Quaternion " << flush;
        //Eigen::Quaterniond _quaternion_e = get_quaternion( system_e );
        // report results
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis("electron", gyro_electron,system_e);
      }
      return;
    }


    void Density2Gyration::AnalyzeGeometry(std::vector<QMAtom*> _atoms){
    
        tools::Elements _elements; 
        double mass=0.0;
        tools::vec centroid = tools::vec(0.0);
        tools::matrix gyration = tools::matrix(0.0);
        std::vector< QMAtom* >::iterator at;
        for (at=_atoms.begin();at<_atoms.end();++at){
            double m = _elements.getMass((*at)->getType());
            const tools::vec & pos =(*at)->getPos();
            mass+= m;
            centroid+=m*pos;
            gyration+=m*(pos|pos);
        }
        centroid/=mass;
        gyration/=mass;
        gyration-=(centroid|centroid);
        Gyrationtensor gyro;
        gyro.mass=mass;
        gyro.centroid=centroid;
        gyro.gyration=gyration;
        tools::matrix::eigensystem_t system;
        gyration.SolveEigensystem(system); 
        //Eigen::Quaterniond _quaternion = get_quaternion( system );
        // report results
        ReportAnalysis( "geometry", gyro,system );   
    }

    void Density2Gyration::ReportAnalysis(string label,Gyrationtensor gyro, tools::matrix::eigensystem_t system){
       
            CTP_LOG(ctp::logINFO, *_log) << "---------------- " << label << " ----------------" << flush;
            CTP_LOG(ctp::logINFO, *_log) << (boost::format("  Norm               = %1$9.4f ") % (gyro.mass) ) << flush;
         
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid x         = %1$9.4f Ang") % (gyro.centroid.getX()*tools::conv::bohr2ang) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid y         = %1$9.4f Ang") % (gyro.centroid.getY()*tools::conv::bohr2ang) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Centroid y         = %1$9.4f Ang") % (gyro.centroid.getZ()*tools::conv::bohr2ang) ) << flush;
            
            double RA2 = tools::conv::bohr2ang  *tools::conv::bohr2ang;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xx = %1$9.4f Ang^2") % (gyro.gyration[0][0]*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xy = %1$9.4f Ang^2") % (gyro.gyration[0][1]*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor xz = %1$9.4f Ang^2") % (gyro.gyration[0][2]*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor yy = %1$9.4f Ang^2") % (gyro.gyration[1][1]*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor yz = %1$9.4f Ang^2") % (gyro.gyration[1][2]*RA2) ) << flush;
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor zz = %1$9.4f Ang^2") % (gyro.gyration[2][2]*RA2) ) << flush;
            
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D1 = %1$9.4f Ang^2") % (system.eigenvalues[0] *RA2) ) << flush;     
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D2 = %1$9.4f Ang^2") % (system.eigenvalues[1]*RA2) ) << flush;     
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Gyration Tensor D3 = %1$9.4f Ang^2") % (system.eigenvalues[2]*RA2) ) << flush;   

            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Radius of Gyration = %1$9.4f Ang") % (std::sqrt(system.eigenvalues[0] +system.eigenvalues[1]+ system.eigenvalues[2])*tools::conv::bohr2ang )) << flush;  
            
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 1 = %1$9.4f ") % (system.eigenvecs[0].getX()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 2 = %1$9.4f ") % (system.eigenvecs[0].getY()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 1 3 = %1$9.4f ") % (system.eigenvecs[0].getZ()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 1 = %1$9.4f ") % (system.eigenvecs[1].getX()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 2 = %1$9.4f ") % (system.eigenvecs[1].getY()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 2 3 = %1$9.4f ") % (system.eigenvecs[1].getZ()) ) << flush; 
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 1 = %1$9.4f ") % (system.eigenvecs[2].getX()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 2 = %1$9.4f ") % (system.eigenvecs[2].getY()) ) << flush;             
            CTP_LOG(ctp::logINFO,*_log) << (boost::format("  Tensor EF Axis 3 3 = %1$9.4f ") % (system.eigenvecs[2].getZ()) ) << flush;             
            return;
    }

    Eigen::Quaterniond Density2Gyration::get_quaternion(const tools::matrix::eigensystem_t& system){
      Eigen::Matrix3d rot = Eigen::Matrix3d::Zero();
      rot<< system.eigenvecs[0].getX(),system.eigenvecs[0].getY(),system.eigenvecs[0].getZ(),
               system.eigenvecs[1].getX(),system.eigenvecs[1].getY(),system.eigenvecs[1].getZ(),
               system.eigenvecs[2].getX(),system.eigenvecs[2].getY(),system.eigenvecs[2].getZ();
       Eigen::Quaterniond quaternion=Eigen::Quaterniond(rot);
        return quaternion;
        
    }
    

}}


