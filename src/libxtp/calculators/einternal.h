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

#ifndef _VOTCA_XTP_EINTERNAL_H
#define _VOTCA_XTP_EINTERNAL_H

#include <votca/xtp/qmcalculator.h>

namespace votca {
namespace xtp {

class EInternal : public QMCalculator {
 public:
  EInternal(){};
  ~EInternal(){};

  std::string Identify() { return "einternal"; }
  void Initialize(tools::Property &options);
  void ParseEnergiesXML(tools::Property &options);
  bool EvaluateFrame(Topology &top);

 private:
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_xX_nN;
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_nX_nN;
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_xN_xX;

  std::map<std::string, QMStateCarrierStorage<bool> > _seg_has_state;

  std::map<std::string, bool> _has_seg;
};

void EInternal::Initialize(tools::Property &options) {

  /* ---- OPTIONS.XML Structure -----
   *
   * <einternal>
   *
   *      <energiesXML>ENERGIES.XML</energiesXML>
   *
   * </einternal>
   *
   */

  this->ParseEnergiesXML(options);
}

void EInternal::ParseEnergiesXML(tools::Property &opt) {

  // update options with the VOTCASHARE defaults
  UpdateWithDefaults(opt, "xtp");
  std::string key = "options." + Identify();

  std::string energiesXML = opt.get(key + ".energiesXML").as<std::string>();

  std::cout << std::endl
            << "... ... Site, reorg. energies from " << energiesXML << ". "
            << std::flush;

  tools::Property alloc;
  tools::load_property_from_xml(alloc, energiesXML);

  /* --- ENERGIES.XML Structure ---
   *
   * <topology>
   *
   *     <molecules>
   *          <molecule>
   *          <name></name>
   *
   *          <segments>
   *
   *              <segment>
   *              <name></name>
   *
   *              <!-- U_sG_sG, s->state, G->geometry !-->
   *
   *              <U_cC_nN_e></U_cC_nN_e>
   *              <U_cC_nN_h></U_cC_nN_h>
   *
   *              <U_nC_nN_e></U_nC_nN_e>
   *              <U_nC_nN_h></U_nC_nN_h>
   *
   *              <U_cN_cC_e></U_cN_cC_e>
   *              <U_cN_cC_h></U_cN_cC_h>
   *
   *              </segment>
   *
   *              <segment>
   *                  ...
   *
   */

  key = "topology.molecules.molecule";
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

      if (segprop->exists("U_cC_nN_e") && segprop->exists("U_nC_nN_e") &&
          segprop->exists("U_cN_cC_e")) {
        QMStateType e(QMStateType::Electron);

        U_xX_nN.setValue(segprop->get("U_cC_nN_e").as<double>(), e);
        U_nX_nN.setValue(segprop->get("U_nC_nN_e").as<double>(), e);
        U_xN_xX.setValue(segprop->get("U_cN_cC_e").as<double>(), e);
        has_state.setValue(true, e);
      }

      if (segprop->exists("U_cC_nN_h") && segprop->exists("U_nC_nN_h") &&
          segprop->exists("U_cN_cC_h")) {

        QMStateType h(QMStateType::Hole);
        U_xX_nN.setValue(segprop->get("U_cC_nN_h").as<double>(), h);
        U_nX_nN.setValue(segprop->get("U_nC_nN_h").as<double>(), h);
        U_xN_xX.setValue(segprop->get("U_cN_cC_h").as<double>(), h);
        has_state.setValue(true, h);
      }

      if (segprop->exists("U_xX_nN_s") && segprop->exists("U_nX_nN_s") &&
          segprop->exists("U_xN_xX_s")) {
        QMStateType s(QMStateType::Singlet);
        U_xX_nN.setValue(segprop->get("U_xX_nN_s").as<double>(), s);
        U_nX_nN.setValue(segprop->get("U_nX_nN_s").as<double>(), s);
        U_xN_xX.setValue(segprop->get("U_xN_xX_s").as<double>(), s);
        has_state.setValue(true, s);
      }
      if (segprop->exists("U_xX_nN_t") && segprop->exists("U_nX_nN_t") &&
          segprop->exists("U_xN_xX_t")) {
        QMStateType t(QMStateType::Triplet);
        U_xX_nN.setValue(segprop->get("U_xX_nN_t").as<double>(), t);
        U_nX_nN.setValue(segprop->get("U_nX_nN_t").as<double>(), t);
        U_xN_xX.setValue(segprop->get("U_xN_xX_t").as<double>(), t);
        has_state.setValue(true, t);
      }
      _seg_has_state[segName] = has_state;
      _seg_U_xX_nN[segName] = U_xX_nN;
      _seg_U_nX_nN[segName] = U_nX_nN;
      _seg_U_xN_xX[segName] = U_xN_xX;

      _has_seg[segName] = has_seg;
    }
  }
}

bool EInternal::EvaluateFrame(Topology &top) {

  int count = 0;
  for (Segment &seg : top.Segments()) {

    std::string segName = seg.getName();

    try {
      _has_seg.at(segName);
    } catch (const std::exception &out_of_range) {
      std::cout << std::endl
                << "... ... WARNING: No energy information for seg [" << segName
                << "]. Skipping... ";
      continue;
    }

    ++count;

    if (_seg_has_state[segName].getValue(QMStateType::Electron)) {
      seg.setU_xX_nN(_seg_U_xX_nN[segName].getValue(QMStateType::Electron),
                     QMStateType::Electron);
      seg.setU_nX_nN(_seg_U_nX_nN[segName].getValue(QMStateType::Electron),
                     QMStateType::Electron);
      seg.setU_xN_xX(_seg_U_xN_xX[segName].getValue(QMStateType::Electron),
                     QMStateType::Electron);
    }

    if (_seg_has_state[segName].getValue(QMStateType::Hole)) {
      seg.setU_xX_nN(_seg_U_xX_nN[segName].getValue(QMStateType::Hole),
                     QMStateType::Hole);
      seg.setU_nX_nN(_seg_U_nX_nN[segName].getValue(QMStateType::Hole),
                     QMStateType::Hole);
      seg.setU_xN_xX(_seg_U_xN_xX[segName].getValue(QMStateType::Hole),
                     QMStateType::Hole);
    }

    if (_seg_has_state[segName].getValue(QMStateType::Singlet)) {
      seg.setU_xX_nN(_seg_U_xX_nN[segName].getValue(QMStateType::Singlet),
                     QMStateType::Singlet);
      seg.setU_nX_nN(_seg_U_nX_nN[segName].getValue(QMStateType::Singlet),
                     QMStateType::Singlet);
      seg.setU_xN_xX(_seg_U_xN_xX[segName].getValue(QMStateType::Singlet),
                     QMStateType::Singlet);
    }

    if (_seg_has_state[segName].getValue(QMStateType::Triplet)) {
      seg.setU_xX_nN(_seg_U_xX_nN[segName].getValue(QMStateType::Triplet),
                     QMStateType::Triplet);
      seg.setU_nX_nN(_seg_U_nX_nN[segName].getValue(QMStateType::Triplet),
                     QMStateType::Triplet);
      seg.setU_xN_xX(_seg_U_xN_xX[segName].getValue(QMStateType::Triplet),
                     QMStateType::Triplet);
    }
  }

  std::cout << std::endl
            << "... ... Read in site, reorg. energies for " << count
            << " segments. " << std::flush;

  return 1;
}

}  // namespace xtp
}  // namespace votca

#endif  //_VOTCA_XTP_EINTERNAL_H
