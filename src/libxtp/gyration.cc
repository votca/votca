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

void Density2Gyration::Initialize( tools::Property& options) {
    string key = Identify();

    std::string statestring=options.ifExistsReturnElseThrowRuntimeError<string>(key + ".state");
    _state.FromString(statestring);
    _dostateonly = options.ifExistsReturnElseReturnDefault<bool>(key + ".difference_to_groundstate",false);
    _gridsize = options.ifExistsReturnElseReturnDefault<string>(key+".gridsize","medium");
    _openmp_threads = options.ifExistsReturnElseReturnDefault<int>(key + ".openmp",1);
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");   
    }


    void Density2Gyration::AnalyzeDensity(Orbitals & orbitals) {
      int threads = 1;
#ifdef _OPENMP
      if (_openmp_threads > 0) omp_set_num_threads(_openmp_threads);
      threads = omp_get_max_threads();
#endif
      CTP_LOG(ctp::logDEBUG, *_log) << "===== Running on " << threads << " threads ===== " << flush;

      std::vector< QMAtom* > Atomlist = orbitals.QMAtoms();
      Eigen::MatrixXd DMAT_tot;
      BasisSet bs;
      bs.LoadBasisSet(orbitals.getDFTbasis());
      AOBasis basis;
      basis.AOBasisFill(bs, Atomlist);
      AnalyzeGeometry(Atomlist);
      std::vector<Eigen::MatrixXd > DMAT;

      // setup numerical integration grid
      NumericalIntegration numway;
      numway.GridSetup(_gridsize, Atomlist, basis);

      if (!_dostateonly) {
        Eigen::MatrixXd DMATGS = orbitals.DensityMatrixFull(_state);
        Gyrationtensor gyro = numway.IntegrateGyrationTensor(DMAT_tot);
        tools::matrix::eigensystem_t system;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro.gyration.SolveEigensystem(system); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating Quaternion " << flush;
        //Eigen::Quaterniond _quaternion = get_quaternion( system );
        // report results
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis(_state.ToLongString(), gyro,system);


      } else{
        // hole density first
        std::vector<Eigen::MatrixXd > DMAT=orbitals.DensityMatrixExcitedState(_state);
        Gyrationtensor gyro_hole = numway.IntegrateGyrationTensor(DMAT[0]);
        tools::matrix::eigensystem_t system_h;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro_hole.gyration.SolveEigensystem(system_h); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis("hole", gyro_hole,system_h);
        
        // electron density
        Gyrationtensor gyro_electron = numway.IntegrateGyrationTensor(DMAT[1]);
        tools::matrix::eigensystem_t system_e;
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Converting to Eigenframe " << flush;
        gyro_electron.gyration.SolveEigensystem(system_e); 
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Reporting " << flush;
        ReportAnalysis("electron", gyro_electron,system_e);
      }
      return;
    }


    void Density2Gyration::AnalyzeGeometry(std::vector<QMAtom*> atoms){
    
        tools::Elements elements; 
        double mass=0.0;
        tools::vec centroid = tools::vec(0.0);
        tools::matrix gyration = tools::matrix(0.0);
        std::vector< QMAtom* >::iterator at;
        for (QMAtom* atom:atoms){
            double m = elements.getMass(atom->getType());
            const tools::vec & pos =atom->getPos();
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

    
    

}}


