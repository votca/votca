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
#ifndef _VOTCA_XTP_ANDERSON__H
#define _VOTCA_XTP_ANDERSON__H

#include <iostream>
#include <memory>
#include <vector>
#include <votca/xtp/eigen.h>
namespace votca {
namespace xtp {

class ANDERSON {
 public:
  Eigen::VectorXd NPAndersonMixing(const double alpha);

  void UpdateOutput(const Eigen::VectorXd &newOutput);
  void UpdateInput(const Eigen::VectorXd &newInput);
  void SetOrder( const Index max_history ) { _max_history = max_history +1; };
  bool Info() { return success; }

 private:
  bool success = true;

  std::vector<Eigen::VectorXd> _input;
  std::vector<Eigen::VectorXd> _output;

  Index _max_history = 25;
  Index _iteration = 0;
};

}  // namespace xtp
}  // namespace votca

#endif
