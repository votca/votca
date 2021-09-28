/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "einternal.h"

namespace votca {
namespace xtp {

void EInternal::ParseOptions(const tools::Property &options) {

  energies_file_ = options.get(".energies_file").as<std::string>();
}

void EInternal::ParseEnergies() {

  std::cout << "\n... ... Site, reorg. energies from " << energies_file_
            << ".\n";

  tools::Property alloc;
  alloc.LoadFromXML(energies_file_);

  std::string key = "topology.molecules.molecule";
  std::vector<tools::Property *> mols = alloc.Select(key);
  for (tools::Property *molprop : mols) {

    key = "segments.segment";
    std::vector<tools::Property *> segs = molprop->Select(key);
    for (tools::Property *segprop : segs) {

      std::string segName = segprop->get("name").as<std::string>();

      bool has_seg = true;

      QMStateCarrierStorage<double> U_xX_nN;
      QMStateCarrierStorage<double> U_nX_nN;
      QMStateCarrierStorage<double> U_xN_xX;
      QMStateCarrierStorage<bool> has_state;
      double eV2hrt = tools::conv::ev2hrt;

      std::vector<QMStateType> types = {QMStateType::Electron,
                                        QMStateType::Hole, QMStateType::Singlet,
                                        QMStateType::Triplet};
      for (QMStateType type : types) {
        std::string u_xX_nN = "U_xX_nN_" + type.ToString();
        std::string u_nX_nN = "U_nX_nN_" + type.ToString();
        std::string u_xN_xX = "U_xN_xX_" + type.ToString();
        if (segprop->exists(u_xX_nN) && segprop->exists(u_nX_nN) &&
            segprop->exists(u_xN_xX)) {
          U_xX_nN.setValue(segprop->get(u_xX_nN).as<double>() * eV2hrt, type);
          U_nX_nN.setValue(segprop->get(u_nX_nN).as<double>() * eV2hrt, type);
          U_xN_xX.setValue(segprop->get(u_xN_xX).as<double>() * eV2hrt, type);
          has_state.setValue(true, type);
        }
      }
      seg_has_state_[segName] = has_state;
      seg_U_xX_nN_[segName] = U_xX_nN;
      seg_U_nX_nN_[segName] = U_nX_nN;
      seg_U_xN_xX_[segName] = U_xN_xX;

      has_seg_[segName] = has_seg;
    }
  }
}

bool EInternal::Evaluate(Topology &top) {

  ParseEnergies();

  Index count = 0;
  for (Segment &seg : top.Segments()) {

    std::string segName = seg.getType();

    if (!has_seg_.count(segName)) {
      std::cout << std::endl
                << "... ... WARNING: No energy information for seg [" << segName
                << "]. Skipping... ";
      continue;
    }

    ++count;

    std::vector<QMStateType> types = {QMStateType::Electron, QMStateType::Hole,
                                      QMStateType::Singlet,
                                      QMStateType::Triplet};
    for (QMStateType type : types) {

      if (seg_has_state_[segName].getValue(type.Type())) {
        seg.setU_xX_nN(seg_U_xX_nN_[segName].getValue(type.Type()),
                       type.Type());
        seg.setU_nX_nN(seg_U_nX_nN_[segName].getValue(type.Type()),
                       type.Type());
        seg.setU_xN_xX(seg_U_xN_xX_[segName].getValue(type.Type()),
                       type.Type());
      }
    }
  }

  std::cout << std::endl
            << "... ... Read in site, reorg. energies for " << count
            << " segments. " << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca
