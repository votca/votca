/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <iomanip>

// Local VOTCA includes
#include "votca/xtp/density_integration.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/staticregion.h"
#include "votca/xtp/aopotential.h"

// Local private VOTCA includes
#include "envcorr.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/eeinteractor.h"


namespace votca {
namespace xtp {

void ENVCORR::ParseOptions(const tools::Property &options) {

  grid_accuracy_ = options.get(".grid").as<std::string>();
  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".input", job_name_ + ".orb");
  outputfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + "_state.dat");
  state_ = options.get(".state").as<QMState>();
  envmode_ = options.get(".envmode").as<std::string>();

  statetype_ = options.get(".statetype").as<std::string>();
  statenumbers_as_string_ =
        options.get(".statenumber").as<std::string>();
}

bool ENVCORR::Run() {

  std::vector<Index> stateindices =
      IndexParser().CreateIndexVector(statenumbers_as_string_);


  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  // get the embedded QM results
  Orbitals orb;
  orb.ReadFromCpt("prodan.orb");
  AOBasis basis = orb.getDftBasis();
  QMMolecule mol = orb.QMAtoms();
/*
  AOOverlap dftOL;
  dftOL.Fill(basis);
  Eigen::MatrixXd OLsqrt = dftOL.Sqrt();
  Eigen::MatrixXd MOs = orb.MOs().eigenvectors();
  Eigen::MatrixXd eps = orb.MOs().eigenvalues().asDiagonal();


  // reconstructing static Hamiltonian
  Eigen::MatrixXd ham_static = OLsqrt.transpose() * OLsqrt *  MOs * eps  * MOs.inverse();
  double eps_homo = (MOs.col(60).transpose() * ham_static * MOs.col(60))(0,0);
  XTP_LOG(Log::error, log_) << TimeStamp()
                                    << " reconstructed static " << eps_homo << std::flush;

  // get MOs from polar KS60
  Orbitals orb_polar;
  orb_polar.ReadFromCpt("prodan_gs.orb");

  Eigen::MatrixXd MOs_polar = orb_polar.MOs().eigenvectors();

  double eps_polar = (MOs_polar.col(60).transpose() * ham_static * MOs_polar.col(60))(0,0);


XTP_LOG(Log::error, log_) << TimeStamp()
                                    << " exp <polar|H_static|polar>  " << eps_polar
<< std::flush;

*/




  bool split_dipoles = true;

  double env_en_GS = 0.0;
  if (envmode_ == "thole"){
      CheckpointFile cpf("groundstate.hdf5", CheckpointAccessLevel::READ);
      CheckpointReader r = cpf.getReader("region_1");
      PolarRegion mmregion(0,log_);
      mmregion.ReadFromCpt(r);
      Index count = 0;

      // get the ground state reference dipoles for Thole
      std::vector<std::unique_ptr<StaticSite>> externalsites_groundstate;
      for ( auto segment : mmregion ){
        for (auto site : segment){
          

          // splitting the dipoles
          if (split_dipoles){
            // Calculate virtual charge positions
            double a = 0.1;  // this is in a0
            double mag_d = site.getDipole().norm();                   // this is in e * a0
            const Eigen::Vector3d dir_d = site.getDipole().normalized();
            const Eigen::Vector3d A = site.getPos() + 0.5 * a * dir_d;
            const Eigen::Vector3d B = site.getPos() - 0.5 * a * dir_d;
            double qA = mag_d / a;
            double qB = -qA;
            if (std::abs(qA) > 1e-12) {
              externalsites_groundstate.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              externalsites_groundstate[count]->setPos(A);
              Vector9d multipole = Vector9d::Zero();  
              multipole(0) = qA;
              externalsites_groundstate[count]->setMultipole(multipole, 0);
              count++;
              externalsites_groundstate.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              externalsites_groundstate[count]->setPos(B);
              multipole(0) = qB;
              externalsites_groundstate[count]->setMultipole(multipole, 0);
              count++;
            }
          } else {
              externalsites_groundstate.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              // copy the induced dipoles to static dipoles, not sure how this works in qmmm though
              Vector9d multipole = Vector9d::Zero();  
              multipole.segment<3>(1) = site.getDipole();
              externalsites_groundstate[count]->setMultipole(multipole, 1);
              count++;
          }          
        }
      }



      // contributions from the nuclear charges
      double E_ext = 0;
      eeInteractor interactor;
      for (const QMAtom& atom : mol) {
        StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
        for (const std::unique_ptr<StaticSite>& site : externalsites_groundstate) {
          if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
            XTP_LOG(Log::error, log_) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << std::flush;
            continue;
          }
          E_ext += interactor.CalcStaticEnergy_site(*site, nucleus);
        }
      }

      //XTP_LOG(Log::error, log_)
      //  << TimeStamp() << " Contribution from nuclear-GSpol " << E_ext
      //  << std::flush;

      // setup the corresponding AOPotential matrix
      AOMultipole dftAOESP_GS;
      dftAOESP_GS.FillPotential(basis, externalsites_groundstate);
      XTP_LOG(Log::error, log_)
        << TimeStamp() << " Filled DFT external multipole potential matrix for the ground state"
        << std::flush;

      // calculate int(rho(r) DeltaV_n(r) d3r)
      QMState state = QMState("n");
      Eigen::MatrixXd dmat = orb.DensityMatrixFull(state) ;

      env_en_GS = dmat.cwiseProduct(dftAOESP_GS.Matrix()).sum();
      XTP_LOG(Log::error, log_)
        << TimeStamp() << " Contribution from DFT (el density      ) in GS env " << env_en_GS  
        << std::flush;

      XTP_LOG(Log::error, log_)
        << TimeStamp() << " Contribution from DFT (             nuc) in GS env " << E_ext  
        << std::flush;

      env_en_GS +=  E_ext;

      XTP_LOG(Log::error, log_)
        << TimeStamp() << " Contribution from DFT (el density + nuc) in GS env " << env_en_GS  
        << std::flush;

for (auto it = stateindices.cbegin() ; it != stateindices.cend(); ++it){
  state = QMState(statetype_ + std::to_string(*it));
  Eigen::MatrixXd dmat_ks = orb.DensityMatrixKSstate(state) ;
  double env_en = dmat_ks.cwiseProduct(dftAOESP_GS.Matrix()).sum();
  XTP_LOG(Log::error, log_) << TimeStamp() << " State " << (*it) << " Single particle energy correction in GS env " << env_en << std::flush;
}
  }

  // get the induced dipoles in the inactive region from MMMM checkpoint file
  for (auto it = stateindices.cbegin() ; it != stateindices.cend(); ++it){
    std::vector<std::unique_ptr<StaticSite>> externalsites;

    if (envmode_ == "thole"){

      CheckpointFile cpf(statetype_ + std::to_string(*it) + ".hdf5", CheckpointAccessLevel::READ);
      CheckpointReader r = cpf.getReader("region_1");
      PolarRegion mmregion(0,log_);
      mmregion.ReadFromCpt(r);
      Index count = 0;
      for ( auto segment : mmregion ){
        for (auto site : segment){
          // splitting the dipoles
          if (split_dipoles){
            // Calculate virtual charge positions
            double a = 0.1;  // this is in a0
            double mag_d = site.getDipole().norm();                   // this is in e * a0
            const Eigen::Vector3d dir_d = site.getDipole().normalized();
            const Eigen::Vector3d A = site.getPos() + 0.5 * a * dir_d;
            const Eigen::Vector3d B = site.getPos() - 0.5 * a * dir_d;
            double qA = mag_d / a;
            double qB = -qA;
            if (std::abs(qA) > 1e-12) {
              externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              externalsites[count]->setPos(A);
              Vector9d multipole = Vector9d::Zero();  
              multipole(0) = qA;
              externalsites[count]->setMultipole(multipole, 0);
              count++;
              externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              externalsites[count]->setPos(B);
              multipole(0) = qB;
              externalsites[count]->setMultipole(multipole, 0);
              count++;
            }
          } else {
              externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
              // copy the induced dipoles to static dipoles, not sure how this works in qmmm though
              Vector9d multipole = Vector9d::Zero();  
              multipole.segment<3>(1) = site.getDipole();
              externalsites[count]->setMultipole(multipole, 1);
              count++;
          }
        }
      }



      XTP_LOG(Log::error, log_)
      << TimeStamp() << " Number of moments " << count
      << std::flush;
    } else if (envmode_ == "dftb"){
      // get the response charges from an mps file
      StaticRegion region(0, log_);
      StaticSegment seg = StaticSegment("", 0);
      std::string mpsfile = statetype_ + "_" + std::to_string(*it) + ".mps";
      seg.LoadFromFile(mpsfile);
      region.push_back(seg);
      for (auto segment : region) {
        for (auto site : segment) {
          externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
        }
      }
    }


    // contributions from the nuclear charges
    double E_ext = 0;
    eeInteractor interactor;
    for (const QMAtom& atom : mol) {
      StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
      for (const std::unique_ptr<StaticSite>& site : externalsites) {
        if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
          XTP_LOG(Log::error, log_) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << std::flush;
           continue;
        }
        E_ext += interactor.CalcStaticEnergy_site(*site, nucleus);
      }
    }

    // setup AOmatrix for dipoles
    AOMultipole dftAOESP;
    dftAOESP.FillPotential(basis, externalsites);

    AOOverlap dftOL;
    dftOL.Fill(basis);
    
    // calculate expectaction value
    // get density matrix of state of interest
   // calculate int(rho(r) DeltaV_n(r) d3r)
      QMState state = QMState("n");
      Eigen::MatrixXd dmat = orb.DensityMatrixFull(state) ;

      env_en_GS = dmat.cwiseProduct(dftAOESP.Matrix()).sum();
      XTP_LOG(Log::error, log_)
        << TimeStamp() <<  " State " << (*it) << " Contribution from DFT (el density      ) in KS env " << env_en_GS  
        << std::flush;

      XTP_LOG(Log::error, log_)
        << TimeStamp() <<  " State " << (*it) << " Contribution from DFT (             nuc) in KS env " << E_ext  
        << std::flush;

      env_en_GS +=  E_ext;

      XTP_LOG(Log::error, log_)
        << TimeStamp() <<  " State " << (*it) << " Contribution from DFT (el density + nuc) in KS env " << env_en_GS  
        << std::flush;

    state = QMState(statetype_ + std::to_string(*it));
    Eigen::MatrixXd dmat_ks = orb.DensityMatrixKSstate(state) ;
    double env_en = dmat_ks.cwiseProduct(dftAOESP.Matrix()).sum();
  
    // try second order corrections
    Eigen::VectorXd  MO = orb.MOs().eigenvectors().col(state.StateIdx());
    double env_en_second = 0.0;
    Eigen::MatrixXd precalc = dftAOESP.Matrix() * MO;

    Index this_state = (*it);
    double QPen = orb.QPpertEnergies()[this_state];
    for ( Index i = 0; i < orb.getGWAmax() ; i++ ) {
      if ( i != this_state ){
        double expval = (orb.MOs().eigenvectors().col(i).transpose() * precalc)(0,0);
        env_en_second += std::pow(expval,2) / ( QPen - orb.QPpertEnergies()(i)) ;
       }
    }



    XTP_LOG(Log::error, log_)
        << TimeStamp() << " State " << (*it) << " Single particle energy correction in state env " << env_en << " and 2nd order " << env_en_second << std::flush;
    }

  return true;

}

}  // namespace xtp
}  // namespace votca
 