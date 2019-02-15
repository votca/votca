/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef VOTCA_XTP_SEGMENTSEXTRACTOR_H
#define VOTCA_XTP_SEGMENTSEXTRACTOR_H

#include <boost/format.hpp>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/propertyiomanipulator.h>

namespace votca {
namespace xtp {

class SegmentsExtractor : public ctp::QMCalculator {
 public:
  SegmentsExtractor(){};
  ~SegmentsExtractor(){};

  std::string Identify() { return "extract.segments"; }
  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
};

void SegmentsExtractor::Initialize(tools::Property *options) { return; }

bool SegmentsExtractor::EvaluateFrame(ctp::Topology *top) {

  // Rigidify std::system (if possible)
  if (!top->Rigidify()) return 0;

  std::string xmlfile = Identify() + ".xml";

  tools::Property state("state", "", "");
  tools::Property &segs = state.add("segments", "");
  tools::Property *next = NULL;

  using boost::format;

  // SEGMENTS
  std::vector<ctp::Segment *>::iterator sit;
  next = &segs;
  for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
    ctp::Segment *seg = *sit;
    tools::Property &segprop = next->add("segment", "");
    segprop.add("id", (format("%1$d") % seg->getId()).str());
    segprop.add("name", seg->getName());
    segprop.add("molecule", seg->getMolecule()->getName());
    segprop.add("xyz",
                (format("%1$+1.4f %2$+1.4f %3$+1.4f") % seg->getPos().getX() %
                 seg->getPos().getY() % seg->getPos().getZ())
                    .str());
    if (seg->hasChrgState(+1)) {
      tools::Property &channel = segprop.add("channel", "");
      channel.setAttribute("type", "hole");
      channel.add("energy_hn",
                  (format("%1$+1.7f") % seg->getSiteEnergy(+1)).str());
      channel.add("lambda_hn",
                  (format("%1$+1.7f") % seg->getU_nC_nN(+1)).str());
      channel.add("lambda_nh",
                  (format("%1$+1.7f") % seg->getU_cN_cC(+1)).str());
      channel.add("occupation_h", (format("%1$+1.7e") % seg->getOcc(+1)).str());
    }
    if (seg->hasChrgState(-1)) {
      tools::Property &channel = segprop.add("channel", "");
      channel.setAttribute("type", "electron");
      channel.add("energy_en",
                  (format("%1$+1.7f") % seg->getSiteEnergy(-1)).str());
      channel.add("lambda_en",
                  (format("%1$+1.7f") % seg->getU_nC_nN(-1)).str());
      channel.add("lambda_ne",
                  (format("%1$+1.7f") % seg->getU_cN_cC(-1)).str());
      channel.add("occupation_e", (format("%1$+1.7e") % seg->getOcc(-1)).str());
    }

    if (tools::globals::verbose) {

      // FRAGMENTS
      tools::Property &fragprop = segprop.add("fragment", "");
      std::vector<ctp::Fragment *>::iterator fit;
      for (fit = seg->Fragments().begin(); fit < seg->Fragments().end();
           ++fit) {
        ctp::Fragment *frag = *fit;
        fragprop.add("id", (format("%1$d") % frag->getId()).str());
        fragprop.add("name", frag->getName());
        fragprop.add("xyz", (format("%1$+1.4f %2$+1.4f %3$+1.4f") %
                             frag->getPos().getX() % frag->getPos().getY() %
                             frag->getPos().getZ())
                                .str());

        // ATOMS
        std::vector<ctp::Atom *>::iterator ait;
        for (ait = frag->Atoms().begin(); ait < frag->Atoms().end(); ++ait) {
          tools::Property &atomprop = fragprop.add("atom", "");
          ctp::Atom *atm = *ait;
          atomprop.add("id", (format("%1$d") % atm->getId()).str());
          atomprop.add("element", atm->getElement());
          atomprop.add("name", atm->getName());
          atomprop.add("weight", (format("%1$1.2f") % atm->getWeight()).str());
          atomprop.add("pos", (format("%1$+1.4f %2$+1.4f %3$+1.4f") %
                               atm->getPos().getX() % atm->getPos().getY() %
                               atm->getPos().getZ())
                                  .str());
          if (atm->HasQMPart()) {
            atomprop.add("qmid", (format("%1$d") % atm->getQMId()).str());
            atomprop.add(
                "qmpos",
                (format("%1$+1.4f %2$+1.4f %3$+1.4f") % atm->getQMPos().getX() %
                 atm->getQMPos().getY() % atm->getQMPos().getZ())
                    .str());
          }
        }
      }
    }
  }

  std::ofstream ofs;
  ofs.open(xmlfile.c_str(), std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + xmlfile);
  }
  ofs << tools::XML << state;
  ofs.close();

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENTSEXTRACTOR_H
