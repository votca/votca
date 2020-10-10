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

class MapChecker : public QMCalculator {
 public:
  MapChecker() = default;

  ~MapChecker() final = default;

  std::string Identify() final { return "mapchecker"; }
  bool WriteToStateFile() const final { return false; }
  void ParseOptions(const tools::Property& user_options) final;
  bool Evaluate(Topology& top) final;

 private:
  std::string AddSteptoFilename(const std::string& filename, Index step) const;
  std::string AddStatetoFilename(const std::string& filename,
                                 QMState state) const;

  std::vector<QMState> StringToStates(const std::string& states_string) const;

  std::string _segmentfile;
  std::string _qmfile;
  std::string _mpfile;
  std::string _mapfile = "";

  std::vector<QMState> _qmstates;
  std::vector<QMState> _mdstates;
};

void MapChecker::ParseOptions(const tools::Property& options) {

  _segmentfile = options.get(".md_pdbfile").as<std::string>();

  _qmfile = options.get(".qm_pdbfile").as<std::string>();

  _mpfile = options.get(".mp_pdbfile").as<std::string>();

  std::string output_qm =
      options.ifExistsReturnElseReturnDefault<std::string>(".qm_states", "");

  _qmstates = StringToStates(output_qm);
  std::string output_md =
      options.ifExistsReturnElseReturnDefault<std::string>(".mp_states", "");
  _mdstates = StringToStates(output_md);
  if (!(_qmstates.empty() && _mdstates.empty())) {
    _mapfile = options.get(".map_file").as<std::string>();
  }
}

std::vector<QMState> MapChecker::StringToStates(
    const std::string& states_string) const {
  std::vector<QMState> result;
  tools::Tokenizer tok_states(states_string, " \t\n");
  std::vector<std::string> states = tok_states.ToVector();
  for (const std::string& s : states) {
    result.push_back(QMState(s));
  }
  return result;
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
  std::string filename = AddSteptoFilename(_segmentfile, top.getStep());
  std::cout << "Writing segments to " << filename << std::endl;
  top.WriteToPdb(filename);

  Logger log;
  log.setReportLevel(Log::current_level);
  log.setCommonPreface("\n... ...");

  QMMapper map(log);

  for (QMState state : _qmstates) {
    map.LoadMappingFile(_mapfile);
    std::string filename_qm = AddStatetoFilename(_qmfile, state);
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
  for (QMState state : _mdstates) {
    mp.LoadMappingFile(_mapfile);
    std::string filename_mp = AddStatetoFilename(_mpfile, state);
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
