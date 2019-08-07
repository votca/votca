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
#ifndef __VOTCA_XTP_NEIGHBORLIST_H
#define __VOTCA_XTP_NEIGHBORLIST_H

#include <boost/format.hpp>
#include <boost/progress.hpp>
#include <votca/tools/globals.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/qmnblist.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {

class Neighborlist : public QMCalculator {
 public:
  std::string Identify() { return "neighborlist"; }
  bool WriteToStateFile() const { return true; }
  void Initialize(tools::Property& options);
  bool EvaluateFrame(Topology& top);

 private:
  int DetClassicalPairs(Topology& top);

  std::vector<std::string> _included_segments;
  std::map<std::string, std::map<std::string, double> > _cutoffs;
  bool _useConstantCutoff;
  double _constantCutoff;
  bool _useExcitonCutoff;
  double _excitonqmCutoff;
};

}  // namespace xtp
}  // namespace votca

#endif
