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
#include "votca/xtp/aopotential.h"

// Local private VOTCA includes
#include "envcorr.h"

namespace votca {
namespace xtp {

void ENVCORR::ParseOptions(const tools::Property &options) {

  grid_accuracy_ = options.get(".grid").as<std::string>();
  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".input", job_name_ + ".orb");
  outputfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + "_state.dat");
  state_ = options.get(".state").as<QMState>();
}

bool ENVCORR::Run() {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  // get the embedded QM results
  Orbitals orb;
  orb.ReadFromCpt(orbfile_);
  AOBasis basis = orb.getDftBasis();
 
  // get the induced dipoles in the inactive region from MMMM checkpoint file
  CheckpointFile cpf("MMMM.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("region_1");
  PolarRegion mmregion(0,log_);
  mmregion.ReadFromCpt(r);
  std::vector<std::unique_ptr<StaticSite>> externalsites;
  Index count = 0;
  for ( auto segment : mmregion ){
    for (auto site : segment){
      externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
      // copy the induced dipoles to static dipoles, not sure how this works in qmmm though
      Vector9d multipole = Vector9d::Zero();  
      multipole.segment<3>(1) = site.getDipole();
      externalsites[count]->setMultipole(multipole, 1);
      count++;
    }
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Number of moments " << count
      << std::flush;

  // setup AOmatrix for dipoles
  AOMultipole dftAOESP;
  dftAOESP.FillPotential(basis, externalsites);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Filled DFT external multipole potential matrix"
      << std::flush;

  // calculate expectaction value
  // get density matrix of state of interest
  QMState state = QMState("ks60");
  Eigen::MatrixXd dmat = orb.DensityMatrixKSstate(state) ;
  double env_en = dmat.cwiseProduct(dftAOESP.Matrix()).sum();

  // try second order corrections
  Eigen::VectorXd  MO = orb.MOs().eigenvectors().col(state.StateIdx());
  double env_en_second = 0.0;
  Eigen::MatrixXd precalc = dftAOESP.Matrix() * MO;
  double QPen = orb.QPpertEnergies()[60];
  for ( Index i = 0; i < orb.getGWAmax() ; i++ ) {
    if ( i != 60 ){
      std::cout << "col " << (orb.MOs().eigenvectors().col(i).transpose() * precalc).cols() << std::endl;
      double expval = (orb.MOs().eigenvectors().col(i).transpose() * precalc)(0,0);
      env_en_second += std::pow(expval,2) / ( QPen - orb.QPpertEnergies()(i)) ;


    }
  }



  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Energy correction " << env_en << " or " << env_en_second
      << std::flush;


  return true;
}

}  // namespace xtp
}  // namespace votca
 