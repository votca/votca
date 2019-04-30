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

#ifndef VOTCA_XTP_OCCUPATIONSEXTRACTOR_H
#define VOTCA_XTP_OCCUPATIONSEXTRACTOR_H

#include <boost/format.hpp>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/propertyiomanipulator.h>

namespace votca {
namespace xtp {

class OccupationsExtractor : public ctp::QMCalculator {
 public:
  OccupationsExtractor(){};
  ~OccupationsExtractor(){};

  std::string Identify() { return "extract.occupations"; }
  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
};

void OccupationsExtractor::Initialize(tools::Property *options) { return; }

bool OccupationsExtractor::EvaluateFrame(ctp::Topology *top) {

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
    segprop.add("xyz",
                (format("%1$+1.4f %2$+1.4f %3$+1.4f") % seg->getPos().getX() %
                 seg->getPos().getY() % seg->getPos().getZ())
                    .str());
    if (seg->hasChrgState(+1)) {
      tools::Property &channel = segprop.add("channel", "");
      channel.setAttribute("type", "hole");
      channel.add("occupation_h", (format("%1$+1.7e") % seg->getOcc(+1)).str());
    }
    if (seg->hasChrgState(-1)) {
      tools::Property &channel = segprop.add("channel", "");
      channel.setAttribute("type", "electron");
      channel.add("occupation_e", (format("%1$+1.7e") % seg->getOcc(-1)).str());
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

#endif  // VOTCA_XTP_OCCUPATIONSEXTRACTOR_H
