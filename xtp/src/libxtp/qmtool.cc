/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#include "votca/xtp/qmtool.h"
#include "votca/xtp/eigen.h"

#include "xtp_libint2.h"

namespace votca {
namespace xtp {
void QMTool::Initialize(const tools::Property& options) {
  job_name_ = options.get("job_name").as<std::string>();
  ParseOptions(options);
}

bool QMTool::Evaluate() {
  libint2::initialize();
  OPENMP::setMaxThreads(nThreads_);
  std::cout << " Using " << OPENMP::getMaxThreads() << " threads" << std::flush;
  bool success = Run();
  libint2::finalize();
  return success;
}

}  // namespace xtp
}  // namespace votca
