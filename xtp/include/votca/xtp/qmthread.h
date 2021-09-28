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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_QMTHREAD_H
#define VOTCA_XTP_QMTHREAD_H

// Standard includes
#include <ctime>
#include <iostream>
#include <string>

// VOTCA includes
#include <votca/tools/globals.h>
#include <votca/tools/thread.h>

// Local VOTCA includes
#include "logger.h"

namespace votca {
namespace xtp {

// ++++++++++++++++++++++++++++++++++++++ //
// Thread class with local string stream //
// ++++++++++++++++++++++++++++++++++++++ //

class QMThread : public tools::Thread {
 public:
  QMThread() = default;
  QMThread(bool maverick) : maverick_(maverick) { ; }
  ~QMThread() override = default;

  Index getId() const { return id_; }
  void setId(Index id) { id_ = id; }
  bool isMaverick() const { return maverick_; }

  Logger& getLogger() { return logger_; }
  void Run(void) override { ; }

 protected:
  Index id_;
  std::stringstream ss_;
  bool maverick_ = false;
  Logger logger_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMTHREAD_H
