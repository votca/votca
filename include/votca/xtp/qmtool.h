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
#ifndef VOTCA_XTP_QMTOOL_H
#define VOTCA_XTP_QMTOOL_H

// Third party includes
#include <boost/format.hpp>

// VOTCA includes

#include <votca/tools/calculator.h>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class QMTool : public tools::Calculator {
 public:
  QMTool() = default;
  ~QMTool() override = default;

  std::string Identify() const override = 0;
  std::string Package() const final { return "xtp"; }
  void Initialize(const tools::Property& options) final;
  bool Evaluate();

 protected:
  virtual bool Run() = 0;
  virtual void ParseOptions(const tools::Property& opt) = 0;
  std::string job_name_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMTOOL_H
