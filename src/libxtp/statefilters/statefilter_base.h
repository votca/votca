/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef __VOTCA_XTP_STATEFILTERBASE_H
#define __VOTCA_XTP_STATEFILTERBASE_H

#include <string>
#include <votca/tools/property.h>
#include <votca/tools/types.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

/**
    \brief Base Class for statefilter


 */

class StateFilter_base {
 public:
  StateFilter_base() = default;
  virtual ~StateFilter_base() = default;

  virtual std::string Identify() const = 0;
  virtual void Initialize(const tools::Property&) = 0;

  virtual void Info(Logger&) const = 0;

  virtual void UpdateHist(const Orbitals&, QMState) = 0;

  virtual std::vector<Index> CalcIndeces(const Orbitals&,
                                         QMStateType) const = 0;

  virtual void WriteToCpt(CheckpointWriter&) = 0;

  virtual void ReadFromCpt(CheckpointReader&) = 0;

 protected:
  std::vector<Index> ReduceAndSortIndecesUp(const Eigen::VectorXd& overlap,
                                            Index offset,
                                            double threshold) const;
  std::vector<Index> ReduceAndSortIndecesDown(const Eigen::VectorXd& overlap,
                                              Index offset,
                                              double threshold) const;

 private:
  template <bool larger>
  std::vector<Index> ReduceAndSortIndeces(const Eigen::VectorXd& overlap,
                                          Index offset, double threshold) const;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_STATEFILTERBASE_H */
