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

#include "erdiabatizationframe.h"

using std::flush;

namespace votca {
namespace xtp {

void ERDiabatizationFrame::ParseOptions(const tools::Property& user_options) {

  log_.setReportLevel(Log::current_level);

  log_.setMultithreading(true);
  log_.setCommonPreface("\n...");

  tools::Property options = user_options;

  orbfile1_ = options.get(".orb_file1").as<std::string>();
  orbfile2_ = options.get(".orb_file2").as<std::string>();

  options_.state_idx_1 = options.get(".state_idx_1").as<Index>();
  options_.state_idx_2 = options.get(".state_idx_2").as<Index>();
  options_.qmtype = options.get(".qmtype").as<std::string>();
  XTP_LOG(Log::error, log_) << "Type : " << options_.qmtype << flush;

  if (options_.state_idx_1 < 1) {
    throw std::runtime_error("State idx 1 must start from 1.");
  } else {
    XTP_LOG(Log::error, log_) << "State 1 : " << options_.state_idx_1 << flush;
  }

  if (options_.state_idx_2 < 1) {
    throw std::runtime_error("State idx 2 must start from 1.");
  } else {
    XTP_LOG(Log::error, log_) << "State 2 : " << options_.state_idx_2 << flush;
  }

  options_.use_RI = options.get(".use_RI").as<bool>();

  XTP_LOG(Log::error, log_) << flush;
}

bool ERDiabatizationFrame::Run() {

  OPENMP::setMaxThreads(nThreads_);

  // set logger
  log_.setReportLevel(Log::current_level);
  // log_.setReportLevel(Log::error);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n...");

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Reading from orbitals from files: " << orbfile1_
      << " and " << orbfile2_ << flush;

  // Get orbitals objects
  Orbitals orbitals1;
  Orbitals orbitals2;

  orbitals1.ReadFromCpt(orbfile1_);
  orbitals2.ReadFromCpt(orbfile2_);

  ERDiabatization ERDiabatization(orbitals1, orbitals2, &log_);

  if (orbitals1.getTDAApprox()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " " << orbfile1_ << "  was done with TDA." << flush;
  }
  if (orbitals2.getTDAApprox()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " " << orbfile2_ << "  was done with TDA. " << flush;
  }

  ERDiabatization.configure(options_);
  ERDiabatization.setUpMatrices();

  // Calculate angle
  double angle = ERDiabatization.Calculate_angle();

  // Calculate the diabatic Hamiltonian
  Eigen::MatrixXd diabatic_H = ERDiabatization.Calculate_diabatic_H(angle);

  // Printing Output
  XTP_LOG(Log::error, log_)
      << "Diabatic Energy 1: " << diabatic_H(0, 0) * votca::tools::conv::hrt2ev
      << flush;
  XTP_LOG(Log::error, log_)
      << "Diabatic Energy 2: " << diabatic_H(1, 1) * votca::tools::conv::hrt2ev
      << flush;
  XTP_LOG(Log::error, log_)
      << "Diabatic Coupling: " << diabatic_H(1, 0) * votca::tools::conv::hrt2ev
      << flush;
  if (std::abs(diabatic_H(1, 0) - diabatic_H(0, 1)) >
      1e-4 * std::abs(diabatic_H(1, 0))) {
    XTP_LOG(Log::error, log_) << "Different offdiagonal "
                              << diabatic_H(0, 1) * votca::tools::conv::hrt2ev
                              << " --- Check carefully!" << flush;
  }
  return true;
}
}  // namespace xtp
}  // namespace votca