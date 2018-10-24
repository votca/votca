/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/espfit.h>
#include <math.h>
#include <votca/tools/constants.h>

namespace votca { namespace xtp {
  using std::flush;
    
void Espfit::Fit2Density(Orbitals& orbitals,const QMState& state, const AOBasis &basis,std::string gridsize) {

    const Eigen::MatrixXd dmat=orbitals.DensityMatrixFull(state);
    // setting up grid
    Grid grid;
    grid.setupCHELPGGrid(orbitals.QMAtoms());
    XTP_LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << grid.getsize() << " points " << flush;

    // Calculating nuclear potential at gridpoints
    AOOverlap overlap;
    overlap.Fill(basis);
    double N_comp=dmat.cwiseProduct(overlap.Matrix()).sum();
   
    NumericalIntegration numway;

    numway.GridSetup(gridsize,orbitals.QMAtoms(),basis);
    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Setup "<<gridsize<<" Numerical Grid with "<<numway.getGridSize()<<" gridpoints."<< flush;
    double N=numway.IntegrateDensity(dmat);
    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Calculated Densities at Numerical Grid, Number of electrons is "<< N << flush;

    if(std::abs(N-N_comp)>0.001){
        XTP_LOG(logDEBUG, *_log) <<"=======================" << flush;
        XTP_LOG(logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid."<< flush;
        N=N_comp;
        XTP_LOG(logDEBUG, *_log) <<"WARNING: Electronnumber set to "<< N << flush;
        XTP_LOG(logDEBUG, *_log) <<"=======================" << flush;
    }

    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    #pragma omp parallel for
    for ( unsigned i = 0 ; i < grid.getsize(); i++){
        grid.getGridValues()(i)=numway.IntegratePotential(grid.getGridPositions()[i]);
    }


    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Electron contribution calculated"  << flush;
    double netcharge=0.0;
    if (!state.isTransition()){
      EvalNuclearPotential(orbitals.QMAtoms(),grid );
      double Znuc=0.0;
       for (const QMAtom& atom:orbitals.QMAtoms()) {
          Znuc += atom.getNuccharge();
        }
      std::cout<<Znuc<<std::endl;
      netcharge=Znuc-N;
    }
    netcharge = std::round(netcharge);
    XTP_LOG(logDEBUG, *_log) << TimeStamp()<< " Netcharge constrained to " << netcharge << flush;

    FitPartialCharges(orbitals,grid, netcharge);
    return;
    }


void Espfit::EvalNuclearPotential(const QMMolecule& atoms, Grid& grid) {
  
    const std::vector< Eigen::Vector3d >& gridpoints = grid.getGridPositions();
    Eigen::VectorXd& gridvalues=grid.getGridValues();
    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points" << flush;

    for (unsigned i = 0; i < gridpoints.size(); i++) {
        for (unsigned j = 0; j < atoms.size(); j++) {
            const Eigen::Vector3d& posatom=atoms[j].getPos();
            double Znuc=atoms[j].getNuccharge();
            double dist_j = (gridpoints[i]-posatom).norm();
            gridvalues(i) += Znuc / dist_j;
        }
    }
    return;
}


void Espfit::Fit2Density_analytic(Orbitals& orbitals,const QMState& state, const AOBasis &basis) {
    // setting up grid
    Grid grid;
    grid.setupCHELPGGrid(orbitals.QMAtoms());
    const Eigen::MatrixXd dmat=orbitals.DensityMatrixFull(state);
    XTP_LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << grid.getsize() << " points " << std::endl;
    // Calculating nuclear potential at gridpoints
    AOOverlap overlap;
    overlap.Fill(basis);
    double N=dmat.cwiseProduct(overlap.Matrix()).sum();

    double netcharge=0.0;
    if (!state.isTransition()){
      EvalNuclearPotential(orbitals.QMAtoms(),grid );
      double Znuc=0.0;
       for (const QMAtom& atom:orbitals.QMAtoms()) {
          Znuc += atom.getNuccharge();
        }
      netcharge=Znuc-N;
    }
    netcharge = std::round(netcharge);

    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    #pragma omp parallel for
    for ( unsigned i = 0 ; i < grid.getsize(); i++){
         AOESP aoesp;
         aoesp.setPosition(grid.getGridPositions()[i]);
         aoesp.Fill(basis);
         grid.getGridValues()(i) -=dmat.cwiseProduct(aoesp.Matrix()).sum();
          }

   FitPartialCharges(orbitals,grid,netcharge);
  
    return;
    }

void Espfit::FitPartialCharges( Orbitals& orbitals,const Grid& grid,double netcharge ){
    const QMMolecule& atomlist=orbitals.QMAtoms();
  const int NoOfConstraints=1+_regionconstraint.size()+_pairconstraint.size();
  const int matrixSize=atomlist.size()+NoOfConstraints;
    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Setting up Matrices for fitting of size "<< matrixSize <<" x " << matrixSize<< flush;

    const std::vector<Eigen::Vector3d >& gridpoints=grid.getGridPositions();
    const Eigen::VectorXd & potential=grid.getGridValues();
    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Using "<< atomlist.size() <<" Fittingcenters and " << gridpoints.size()<< " Gridpoints."<< flush;

    Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(matrixSize,matrixSize);
    Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(matrixSize);
    // setting up _Amat
    #pragma omp parallel for
    for ( unsigned i =0 ; i < atomlist.size(); i++){
        for ( unsigned j=i; j<atomlist.size(); j++){
            for ( unsigned k=0; k < gridpoints.size(); k++){
                double dist_i = (atomlist[i].getPos()-gridpoints[k]).norm();
                double dist_j = (atomlist[j].getPos()-gridpoints[k]).norm();

                 Amat(i,j) += 1.0/dist_i/dist_j;
            }
            Amat(j,i) = Amat(i,j);
        }
    }

     // setting up Bvec
    #pragma omp parallel for
    for ( unsigned i =0 ; i < atomlist.size(); i++){
        for ( unsigned k=0; k < gridpoints.size(); k++){
                double dist_i = (atomlist[i].getPos()-gridpoints[k]).norm();
                Bvec(i) += potential(k)/dist_i;
        }
       }
    //Total charge constraint
    for ( unsigned i =0 ; i < atomlist.size()+1; i++){
      Amat(i,atomlist.size()) = 1.0;
      Amat(atomlist.size(),i) = 1.0;
    }
    Amat(atomlist.size(),atomlist.size()) = 0.0;
    Bvec(atomlist.size()) = netcharge; //netcharge!!!!
     
     
    // Pairconstraint
     for (unsigned i=0;i<_pairconstraint.size();i++){
         const std::pair<int,int>& pair=_pairconstraint[i];
         Amat(pair.first,atomlist.size()+1+i)=1.0;
         Amat(atomlist.size()+1+i,pair.first)=1.0;
         Amat(pair.second,atomlist.size()+1+i)=-1.0;
         Amat(atomlist.size()+1+i,pair.second)=-1.0;
     }
     
     //Regionconstraint
     for (unsigned i=0;i<_regionconstraint.size();i++){
         const ConstraintRegion& reg=_regionconstraint[i];
         for (const int& index:reg.atomindices){
             Amat(index,atomlist.size()+i+1+_pairconstraint.size())=1.0;
             Amat(atomlist.size()+i+1+_pairconstraint.size(),index)=1.0;
         }
         Bvec(atomlist.size()+i+1+_pairconstraint.size())=reg.charge;
     }

    XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Solving linear Equation "<< flush;
    Eigen::VectorXd charges;
    if(_do_svd){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd;
      svd.setThreshold(_conditionnumber);
      svd.compute(Amat,Eigen::ComputeThinU | Eigen::ComputeThinV);
      charges=svd.solve(Bvec);
      XTP_LOG(logDEBUG, *_log) << TimeStamp() << " SVD Done. "<<flush;
      if((Bvec.size()-svd.nonzeroSingularValues())!=0){
        XTP_LOG(logDEBUG, *_log) << TimeStamp() <<Bvec.size()-svd.nonzeroSingularValues()<<" Sites could not be fitted and are set to zero."<< flush;
      }
    }
    else{
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(Amat);
      if(QR.info()!=Eigen::ComputationInfo::Success){
        throw std::runtime_error("Espfit: Solving the constrained equation failed. Maybe try SVD.");
      }
      charges=QR.solve(Bvec);
      XTP_LOG(logDEBUG, *_log) << TimeStamp() << " Solved linear least square fit ."<< flush;
    }
    //remove constraints from charges
    charges.conservativeResize(atomlist.size());
    PolarSegment seg=PolarSegment(orbitals.QMAtoms().getName(),orbitals.QMAtoms().getId());

    XTP_LOG(logDEBUG, *_log) << " Sum of fitted charges: " << charges.sum() << flush;
    for (int i=0;i<atomlist.size();i++){
        seg.push_back(PolarSite(atomlist[i],charges(i)));
    }
    orbitals.Multipoles()=seg;
    // get RMSE
    double rmse = 0.0;
    double totalPotSq = 0.0;
    for ( unsigned k=0 ; k < gridpoints.size(); k++ ){
        double temp = 0.0;
        for ( const PolarSite& atom:orbitals.Multipoles()){
            double dist = (gridpoints[k]-atom.getPos()).norm();
            temp += atom.getCharge()/dist;
        }
        rmse += (potential(k) - temp)*(potential(k) - temp);
        totalPotSq += potential(k)*potential(k);
    }
    XTP_LOG(logDEBUG, *_log) << " RMSE of fit:  " << sqrt(rmse/gridpoints.size()) << flush;
    XTP_LOG(logDEBUG, *_log) << " RRMSE of fit: " << sqrt(rmse/totalPotSq) << flush;

    return;
   }


}}
