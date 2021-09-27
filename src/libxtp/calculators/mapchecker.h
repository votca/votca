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
#ifndef VOTCA_XTP_MAPCHECKER_H
#define VOTCA_XTP_MAPCHECKER_H

// VOTCA includes
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/segmentmapper.h"

namespace votca {
namespace xtp {

class MapChecker final : public QMCalculator {
 public:
  MapChecker() = default;

  ~MapChecker() = default;

  std::string Identify() const { return "mapchecker"; }
  bool WriteToStateFile() const { return false; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  std::string AddSteptoFilename(const std::string& filename, Index step) const;
  std::string AddStatetoFilename(const std::string& filename,
                                 QMState state) const;

  std::vector<QMState> StringToStates(const std::string& states_string) const;

  std::string segmentfile_;
  std::string qmfile_;
  std::string mpfile_;
  std::string mapfile_ = "";

  std::vector<QMState> qmstates_;
  std::vector<QMState> mdstates_;
};

void MapChecker::ParseOptions(const tools::Property& options) {

  segmentfile_ = options.get(".md_pdbfile").as<std::string>();

  qmfile_ = options.get(".qm_pdbfile").as<std::string>();

  mpfile_ = options.get(".mp_pdbfile").as<std::string>();

  qmstates_ = options.get(".qm_states").as<std::vector<QMState>>();

  mdstates_ = options.get(".mp_states").as<std::vector<QMState>>();
  if (!(qmstates_.empty() && mdstates_.empty())) {
    mapfile_ = options.get(".map_file").as<std::string>();
  }
}

std::string MapChecker::AddStatetoFilename(const std::string& filename,
                                           QMState state) const {
  std::string base = tools::filesystem::GetFileBase(filename);
  std::string fileending = tools::filesystem::GetFileExtension(filename);
  std::string filename_comp = base + "_" + state.ToString() + "." + fileending;
  return filename_comp;
}

bool MapChecker::Evaluate(Topology& top) {
  std::cout << std::endl;
  std::string filename = AddSteptoFilename(segmentfile_, top.getStep());
  std::cout << "Writing segments to " << filename << std::endl;
  top.WriteToPdb(filename);

  Logger log;
  log.setReportLevel(Log::current_level);
  log.setCommonPreface("\n... ...");

  QMMapper map(log);

  for (QMState state : qmstates_) {
    map.LoadMappingFile(mapfile_);
    std::string filename_qm = AddStatetoFilename(qmfile_, state);
    csg::PDBWriter qmwriter;
    std::string filename_qm_state =
        AddSteptoFilename(filename_qm, top.getStep());
    std::cout << "Writing qmmolecules to " << filename_qm_state << std::endl;
    qmwriter.Open(filename_qm_state, false);
    qmwriter.WriteHeader("QMFrame:" + std::to_string(top.getStep()) +
                         " state:" + state.ToString());
    qmwriter.WriteBox(top.getBox() * tools::conv::bohr2ang);
    for (const Segment& seg : top.Segments()) {
      QMMolecule mol = map.map(seg, state);
      qmwriter.WriteContainer(mol);
    }
    qmwriter.Close();
  }

  PolarMapper mp(log);
  for (QMState state : mdstates_) {
    mp.LoadMappingFile(mapfile_);
    std::string filename_mp = AddStatetoFilename(mpfile_, state);
    csg::PDBWriter mpwriter;
    std::string filename_mp_state =
        AddSteptoFilename(filename_mp, top.getStep());
    std::cout << "Writing polarsegments to " << filename_mp_state << std::endl;
    mpwriter.Open(filename_mp_state, false);
    mpwriter.WriteHeader("MPFrame:" + std::to_string(top.getStep()) +
                         " state:" + state.ToString());
    mpwriter.WriteBox(top.getBox() * tools::conv::bohr2ang);
    for (const Segment& seg : top.Segments()) {
      PolarSegment mol = mp.map(seg, state);
      mpwriter.WriteContainer(mol);
    }
    mpwriter.Close();
  }

  return true;
}

std::string MapChecker::AddSteptoFilename(const std::string& filename,
                                          Index step) const {
  std::string base = tools::filesystem::GetFileBase(filename);
  std::string fileending = tools::filesystem::GetFileExtension(filename);
  std::string filename_comp =
      base + "_step_" + std::to_string(step) + "." + fileending;
  return filename_comp;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MAPCHECKER_H
