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
#ifndef VOTCA_XTP_STATE_SAVER_H
#define VOTCA_XTP_STATE_SAVER_H

#include <map>
#include <stdio.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {

class StateSaver {
 public:
  StateSaver(std::string file) : _hdf5file(file){};

  void WriteFrame(const Topology &top);

  Topology ReadFrame(int frameid) const;

  std::vector<int> getFrames() const;

 private:
  bool TopStepisinFrames(int frameid) const;

  std::string _hdf5file;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATE_SAVER_H
