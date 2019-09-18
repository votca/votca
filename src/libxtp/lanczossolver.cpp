/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <iostream>
#include <stdexcept>

#include <votca/xtp/lanczossolver.h>
#include <votca/xtp/eigen.h>
#include <Spectra/SymEigsSolver.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

using namespace std;

LanczosSolver::LanczosSolver(Logger &log) : _log(log) {}


void LanczosSolver::PrintTiming(
    const std::chrono::time_point<std::chrono::system_clock> &start) const {
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << "-----------------------------------" << flush;
  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << "- Lanczos ran for "
                               << elapsed_time.count() << "secs." << flush;
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << "-----------------------------------" << flush;
}

void LanczosSolver::PrintOptions(int op_size) const {

  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Lanczos Solver using " << OPENMP::getMaxThreads()
      << " threads." << flush;

  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Matrix size : " << op_size << 'x' << op_size << flush;
}


}  // namespace xtp
}  // namespace votca