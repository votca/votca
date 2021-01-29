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

  std::string key = "erdiabatization";

  _log.setReportLevel(Log::current_level);

  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  tools::Property options = user_options;
     
  _orbfile = options.get(".orb_file").as<std::string>();


  _options.state_idx_1 = options.get(".state_idx_1").as<Index>();
  _options.state_idx_2 = options.get(".state_idx_2").as<Index>();
  std::string qmtype = options.get(".qmtype").as<std::string>();
  _qmtype.FromString(qmtype);

  XTP_LOG(Log::error, _log) << "Type : " << qmtype  << flush;
  XTP_LOG(Log::error, _log) << "State 1 : " <<  _options.state_idx_1 << flush;
  XTP_LOG(Log::error, _log) << "State 2 : " <<  _options.state_idx_2 << flush;

  XTP_LOG(Log::error, _log) << flush;

};

bool ERDiabatizationFrame::Run() {

  OPENMP::setMaxThreads(_nThreads);

  // set logger

  _log.setReportLevel(Log::error);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Reading from orbitals from file: " << _orbfile
      << flush;

  // Get orbitals object
  Orbitals orbitals;

  orbitals.ReadFromCpt(_orbfile);

  ERDiabatization ERDiabatization(orbitals, &_log);

  ERDiabatization.configure(_options);

  ERDiabatization.setUpMatrices();

  
  XTP_LOG(Log::error, _log) << TimeStamp() << " Started ER Diabatization " << flush;
  
  Eigen::VectorXd results = ERDiabatization.CalculateER(orbitals,_qmtype);

  XTP_LOG(Log::error, _log) << TimeStamp() << " Calculation done " << flush;

        
  for (Index n = 1; n < 101; n++) {
    std::cout << n << "  " << results(n-1) << std::endl;
  }
  return true;
}
}  // namespace xtp
}  // namespace votca