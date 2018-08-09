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
#include <votca/ctp/logger.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/aomatrix.h>

//#include <boost/progress.hpp>

#include <math.h>
#include <votca/tools/constants.h>



namespace votca { namespace xtp {
  using std::flush;
    
void Espfit::Fit2Density(std::vector< QMAtom* >& atomlist,const Eigen::MatrixXd &dmat,const AOBasis &basis,std::string gridsize) {


    // setting up grid
    Grid grid;
    grid.setupCHELPGGrid(atomlist);
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<  " Done setting up CHELPG grid with " << grid.getsize() << " points " << flush;

    // Calculating nuclear potential at gridpoints
    AOOverlap overlap;
    overlap.Fill(basis);
    double N_comp=dmat.cwiseProduct(overlap.Matrix()).sum();
   
    NumericalIntegration numway;

    numway.GridSetup(gridsize,atomlist,basis);
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Setup "<<gridsize<<" Numerical Grid with "<<numway.getGridSize()<<" gridpoints."<< flush;
    double N=numway.IntegrateDensity(dmat);
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculated Densities at Numerical Grid, Number of electrons is "<< N << flush;

    if(std::abs(N-N_comp)>0.001){
        CTP_LOG(ctp::logDEBUG, *_log) <<"=======================" << flush;
        CTP_LOG(ctp::logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid."<< flush;
        N=N_comp;
        CTP_LOG(ctp::logDEBUG, *_log) <<"WARNING: Electronnumber set to "<< N << flush;
        CTP_LOG(ctp::logDEBUG, *_log) <<"=======================" << flush;
    }

    double netcharge=getNetcharge( atomlist,N );

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    //boost::progress_display show_progress( _grid.getsize() );
    #pragma omp parallel for
    for ( unsigned i = 0 ; i < grid.getsize(); i++){
        grid.getGridValues()(i)=numway.IntegratePotential(grid.getGridPositions()[i]);
    }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Electron contribution calculated"  << flush;
    if (!_do_Transition){
      EvalNuclearPotential(  atomlist,  grid );
    }

    FitPartialCharges(atomlist,grid, netcharge);  
    return;
    }


void Espfit::EvalNuclearPotential(const std::vector< QMAtom* >& atoms, Grid& grid) {
  
    const std::vector< tools::vec >& _gridpoints = grid.getGridPositions();
    Eigen::VectorXd& _gridvalues=grid.getGridValues();
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points" << flush;

    for (unsigned i = 0; i < _gridpoints.size(); i++) {
        for (unsigned j = 0; j < atoms.size(); j++) {
            const tools::vec& posatom=atoms[j]->getPos();
            double Znuc=atoms[j]->getNuccharge();
            double dist_j = tools::abs(_gridpoints[i]-posatom);
            _gridvalues(i) += Znuc / dist_j;
        }
    }
    return;
}

double Espfit::getNetcharge(const std::vector< QMAtom* >& atoms, double N ){
    double netcharge = 0.0;
      if (std::abs(N) < 0.05) {
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp()<< " Number of Electrons is " << N << " transitiondensity is used for fit" << flush;
        _do_Transition = true;
      } else {
        double Znuc = 0.0;
        for (unsigned j = 0; j < atoms.size(); j++) {
          Znuc += atoms[j]->getNuccharge();
        }

        if (std::abs(Znuc - N) < 4) {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp()<< " Number of Electrons minus Nucleus charge is " << Znuc - N << "." << flush;
        } else {
          CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp()<< " Warning: Your molecule is highly ionized." << flush;
        }
        netcharge = Znuc - N;
        _do_Transition = false;
      }

      netcharge = std::round(netcharge);
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp()<< " Netcharge constrained to " << netcharge << flush;

      return netcharge;
}



void Espfit::Fit2Density_analytic(std::vector< QMAtom* >& atomlist,const Eigen::MatrixXd &dmat,const AOBasis &basis) {
    // setting up grid
    Grid grid;
    grid.setupCHELPGGrid(atomlist);

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<  " Done setting up CHELPG grid with " << grid.getsize() << " points " << std::endl;
    // Calculating nuclear potential at gridpoints
    AOOverlap overlap;
    overlap.Fill(basis);
    double N_comp=dmat.cwiseProduct(overlap.Matrix()).sum();

    double netcharge=getNetcharge( atomlist,N_comp );
    if(!_do_Transition){
        EvalNuclearPotential(atomlist, grid);
    }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    #pragma omp parallel for
    for ( unsigned i = 0 ; i < grid.getsize(); i++){
         AOESP aoesp;
         aoesp.setPosition(grid.getGridPositions()[i]);
         aoesp.Fill(basis);
         grid.getGridValues()(i) -=dmat.cwiseProduct(aoesp.Matrix()).sum();
          }

   FitPartialCharges(atomlist,grid,netcharge);
  
    return;
    }

void Espfit::FitPartialCharges( std::vector< QMAtom* >& atomlist,const Grid& grid,double netcharge ){
  
  const int NoOfConstraints=1+_regionconstraint.size()+_pairconstraint.size();
  const int matrixSize=atomlist.size()+NoOfConstraints;
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Setting up Matrices for fitting of size "<< matrixSize <<" x " << matrixSize<< flush;

    const std::vector< tools::vec >& gridpoints=grid.getGridPositions();
    const Eigen::VectorXd & potential=grid.getGridValues();
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Using "<< atomlist.size() <<" Fittingcenters and " << gridpoints.size()<< " Gridpoints."<< flush;

    Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(matrixSize,matrixSize);
    Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(matrixSize);
    // setting up _Amat
    #pragma omp parallel for
    for ( unsigned _i =0 ; _i < atomlist.size(); _i++){
        for ( unsigned _j=_i; _j<atomlist.size(); _j++){
            for ( unsigned _k=0; _k < gridpoints.size(); _k++){
                double dist_i = tools::abs(atomlist[_i]->getPos()-gridpoints[_k]);
                double dist_j = tools::abs(atomlist[_j]->getPos()-gridpoints[_k]);

                 Amat(_i,_j) += 1.0/dist_i/dist_j;
            }
            Amat(_j,_i) = Amat(_i,_j);
        }
    }

     // setting up Bvec
    #pragma omp parallel for
    for ( unsigned _i =0 ; _i < atomlist.size(); _i++){
        for ( unsigned _k=0; _k < gridpoints.size(); _k++){
                double dist_i = tools::abs(atomlist[_i]->getPos()-gridpoints[_k]);
                Bvec(_i) += potential(_k)/dist_i;
        }
       }
    //Total charge constraint
    for ( unsigned _i =0 ; _i < atomlist.size()+1; _i++){
      Amat(_i,atomlist.size()) = 1.0;
      Amat(atomlist.size(),_i) = 1.0;
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
         const region& reg=_regionconstraint[i];
         for (const int& index:reg.atomindices){
             Amat(index,atomlist.size()+i+1+_pairconstraint.size())=1.0;
             Amat(atomlist.size()+i+1+_pairconstraint.size(),index)=1.0;
         }
         Bvec(atomlist.size()+i+1+_pairconstraint.size())=reg.charge;
     }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Solving linear Equation "<< flush;
    Eigen::VectorXd charges;
    if(_do_svd){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd;
      svd.setThreshold(_conditionnumber);
      svd.compute(Amat,Eigen::ComputeThinU | Eigen::ComputeThinV);
      charges=svd.solve(Bvec);
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " SVD Done. "<<Bvec.size()-svd.nonzeroSingularValues()<<" Sites could not be fitted and are set to zero."<< flush;
    }
    else{
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(Amat);
      if(QR.info()!=Eigen::ComputationInfo::Success){
        throw std::runtime_error("Espfit: Solving the constrained equation failed. Maybe try SVD.");
      }
      charges=QR.solve(Bvec);
      CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Solved linear least square fit ."<< flush;
    }
    //remove constraint
    charges.conservativeResize(atomlist.size());
   
    CTP_LOG(ctp::logDEBUG, *_log) << " Sum of fitted charges: " << charges.sum() << flush;
    for (unsigned i=0;i<atomlist.size();i++){
      atomlist[i]->setPartialcharge(charges(i));
    }
    
    // get RMSE
    double rmse = 0.0;
    double totalPotSq = 0.0;
    for ( unsigned _k=0 ; _k < gridpoints.size(); _k++ ){
        double temp = 0.0;
        for ( const QMAtom* atom:atomlist){
            double dist =  tools::abs(gridpoints[_k]-atom->getPos());
            temp += atom->getPartialcharge()/dist;
        }
        rmse += (potential(_k) - temp)*(potential(_k) - temp);
        totalPotSq += potential(_k)*potential(_k);
    }
    CTP_LOG(ctp::logDEBUG, *_log) << " RMSE of fit:  " << sqrt(rmse/gridpoints.size()) << flush;
    CTP_LOG(ctp::logDEBUG, *_log) << " RRMSE of fit: " << sqrt(rmse/totalPotSq) << flush;

    return;
   }


}}
