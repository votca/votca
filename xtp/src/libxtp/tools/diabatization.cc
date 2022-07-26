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

#include "diabatization.h"

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void Diabatization::ParseOptions(const tools::Property& user_options) {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n...");

  tools::Property options = user_options;

  // getting diabatization method and checking validity
  method_ = options.get(".method").as<std::string>();
  if (method_ == "er") {
    XTP_LOG(Log::error, log_) << "Method : Edminston-Rudenberg" << flush;
  } else if (method_ == "gmh") {
    XTP_LOG(Log::error, log_) << "Method : Generalized Mulliken-Hush" << flush;
  } else if (method_ == "fcd") {
    XTP_LOG(Log::error, log_) << "Method : Fragment Charge Difference" << flush;
  } else {
    throw std::runtime_error("Diabatization method unknown!");
  }

  // getting orbfiles
  orbfile1_ = options.get(".orb_file1").as<std::string>();
  orbfile2_ = options.get(".orb_file2").as<std::string>();

  // getting state options and checking validity
  state_idx_1_ = options.get(".state_idx_1").as<Index>();
  state_idx_2_ = options.get(".state_idx_2").as<Index>();
  qmtype_ = options.get(".qmtype").as<std::string>();
  XTP_LOG(Log::error, log_) << "Type : " << qmtype_ << flush;

  if (state_idx_1_ < 1) {
    throw std::runtime_error("State idx 1 must start from 1.");
  } else {
    XTP_LOG(Log::error, log_) << "State 1 : " << state_idx_1_ << flush;
  }

  if (state_idx_2_ < 1) {
    throw std::runtime_error("State idx 2 must start from 1.");
  } else {
    XTP_LOG(Log::error, log_) << "State 2 : " << state_idx_2_ << flush;
  }

  useRI_ = options.get(".use_RI").as<bool>();

  if (options.exists(".fragments")) {
    std::vector<tools::Property*> prop_region =
        options.Select("fragments.fragment");
    Index index = 0;
    for (tools::Property* prop : prop_region) {
      std::string indices = prop->get("indices").as<std::string>();
      fragments_.push_back(QMFragment<BSE_Population>(index, indices));
      index++;
    }
  }

  XTP_LOG(Log::error, log_) << flush;
}

bool Diabatization::Run() {

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

  if (orbitals1.getTDAApprox()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " " << orbfile1_ << "  was done with TDA." << flush;
  }
  if (orbitals2.getTDAApprox()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " " << orbfile2_ << "  was done with TDA. " << flush;
  }

  if (method_ == "er") {
    ERDiabatization ERDiabatization(orbitals1, orbitals2, &log_, state_idx_1_,
                                    state_idx_2_, qmtype_, useRI_);
    ERDiabatization.configure();
    ERDiabatization.setUpMatrices();

    // Calculate optimal mixing angle
    double angle = ERDiabatization.Calculate_angle();

    // Calculate the diabatic Hamiltonian
    Eigen::MatrixXd diabatic_H = ERDiabatization.Calculate_diabatic_H(angle);

    // Printing Output
    XTP_LOG(Log::error, log_)
        << format("Diabatic Energy 1: %1$+1.12f eV") %
               (diabatic_H(0, 0) * votca::tools::conv::hrt2ev)
        << flush;
    XTP_LOG(Log::error, log_)
        << format("Diabatic Energy 2: %1$+1.12f eV") %
               (diabatic_H(1, 1) * votca::tools::conv::hrt2ev)
        << flush;
    XTP_LOG(Log::error, log_)
        << format("Diabatic Coupling: %1$+1.12f eV") %
               (diabatic_H(1, 0) * votca::tools::conv::hrt2ev)
        << flush;
    if (std::abs(diabatic_H(1, 0) - diabatic_H(0, 1)) >
        1e-4 * std::abs(diabatic_H(1, 0))) {
      XTP_LOG(Log::error, log_) << "Different offdiagonal "
                                << diabatic_H(0, 1) * votca::tools::conv::hrt2ev
                                << " --- Check carefully!" << flush;
    }
  } else if (method_ == "gmh") {

    GMHDiabatization GMHDiabatization(orbitals1, orbitals2, &log_, state_idx_1_,
                                      state_idx_2_, qmtype_);

    GMHDiabatization.configure();

    std::pair<double, double> coupling = GMHDiabatization.calculate_coupling();
    XTP_LOG(Log::error, log_)
        << format("Diabatic Coupling: %1$+1.12f eV") %
               (coupling.first * votca::tools::conv::hrt2ev)
        << flush;
    XTP_LOG(Log::error, log_)
        << format("Diabatic Coupling with CT axis projection: %1$+1.12f eV") %
               (coupling.second * votca::tools::conv::hrt2ev)
        << flush;
  } else if (method_ == "fcd") {

    // check if fragments are empty
    if (fragments_.size() == 0) {
      throw std::runtime_error("Fragments are undefined in FCD!");
    }

    FCDDiabatization FCDDiabatization(orbitals1, orbitals2, &log_, state_idx_1_,
                                      state_idx_2_, qmtype_, fragments_);

    FCDDiabatization.configure();

    double coupling = FCDDiabatization.calculate_coupling();
    XTP_LOG(Log::error, log_) << format("Diabatic Coupling: %1$+1.12f eV") %
                                     (coupling * votca::tools::conv::hrt2ev)
                              << flush;
  }
  return true;
}
}  // namespace xtp
}  // namespace votca