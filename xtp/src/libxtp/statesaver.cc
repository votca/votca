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

// Third party includes
#include <boost/interprocess/sync/file_lock.hpp>

// VOTCA includes
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/statesaver.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {

std::vector<Index> StateSaver::getFrames() const {
  CheckpointFile cpf(hdf5file_, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader();
  std::vector<Index> frames = std::vector<Index>{};
  try {
    r(frames, "frames");
  } catch (std::runtime_error&) {
    ;
  }
  return frames;
}
void StateSaver::WriteFrame(const Topology& top) {
  if (!tools::filesystem::FileExists(hdf5file_)) {
    std::cout << "Creating statefile " << hdf5file_ << std::endl;
    CheckpointFile cpf(hdf5file_, CheckpointAccessLevel::CREATE);
  }
  boost::interprocess::file_lock flock(hdf5file_.c_str());
  flock.lock();
  if (!TopStepisinFrames(top.getStep())) {
    std::vector<Index> frames = this->getFrames();
    frames.push_back(top.getStep());
    CheckpointFile cpf(hdf5file_, CheckpointAccessLevel::MODIFY);
    CheckpointWriter w = cpf.getWriter();
    w(frames, "frames");
    std::cout << "Frame with id " << top.getStep() << " was not in statefile "
              << hdf5file_ << " ,adding it now." << std::endl;
  }

  CheckpointFile cpf(hdf5file_, CheckpointAccessLevel::MODIFY);
  CheckpointWriter w = cpf.getWriter("/frame_" + std::to_string(top.getStep()));
  top.WriteToCpt(w);

  flock.unlock();

  std::cout << "Wrote MD topology (step = " << top.getStep()
            << ", time = " << top.getTime() << ") to " << hdf5file_
            << std::endl;
  std::cout << "... ";

  std::cout << ". " << std::endl;
  return;
}

Topology StateSaver::ReadFrame(Index frameid) const {
  if (!tools::filesystem::FileExists(hdf5file_)) {
    throw std::runtime_error("Statefile " + hdf5file_ + " does not exist.");
  }
  std::cout << "Import MD Topology (i.e. frame " << frameid << ")"
            << " from " << hdf5file_ << std::endl;
  std::cout << "...";
  boost::interprocess::file_lock flock(hdf5file_.c_str());
  flock.lock();
  if (!TopStepisinFrames(frameid)) {
    throw std::runtime_error("Frame with id " + std::to_string(frameid) +
                             " is not in statefile.");
  }

  CheckpointFile cpf(hdf5file_, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("/frame_" + std::to_string(frameid));
  Topology top;
  top.ReadFromCpt(r);
  flock.unlock();
  std::cout << ". " << std::endl;
  return top;
}

bool StateSaver::TopStepisinFrames(Index frameid) const {
  std::vector<Index> frames = this->getFrames();
  return std::find(frames.begin(), frames.end(), frameid) != frames.end();
}

}  // namespace xtp
}  // namespace votca
