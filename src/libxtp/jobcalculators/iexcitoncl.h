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

#pragma once
#ifndef VOTCA_XTP_IEXCITONCL_H
#define VOTCA_XTP_IEXCITONCL_H

// Third party includes
#include <boost/filesystem.hpp>
#include <sys/stat.h>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/parallelxjobcalc.h"
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

/**
 * \brief Evaluates Transition Charge distributions classically
 *
 * Evaluates the electrostatic classical coupling between molecules in
 * their excited states.
 * Callname: iexcitoncl
 */

class IEXCITON final : public ParallelXJobCalc<std::vector<Job> > {
 public:
  std::string Identify() const { return "iexcitoncl"; }

  Job::JobResult EvalJob(const Topology &top, Job &job, QMThread &opThread);

  void WriteJobFile(const Topology &top);
  void ReadJobFile(Topology &top);

 protected:
  void ParseSpecificOptions(const tools::Property &user_options);

 private:
  QMState GetElementFromMap(const std::string &elementname) const;
  std::map<std::string, QMState> FillParseMaps(const std::string &Mapstring);
  double cutoff_;
  std::map<std::string, QMState> statemap_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IEXCITONCL_H
