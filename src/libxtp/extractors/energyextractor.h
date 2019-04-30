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

#ifndef VOTCA_XTP_ENERGYEXTRACTOR_H
#define VOTCA_XTP_ENERGYEXTRACTOR_H

#include <boost/format.hpp>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/propertyiomanipulator.h>

namespace votca {
namespace xtp {

class EnergyExtractor : public ctp::QMCalculator {
 public:
  EnergyExtractor(){};
  ~EnergyExtractor(){};

  std::string Identify() { return "extract.energy"; }
  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
};

void EnergyExtractor::Initialize(tools::Property *options) { return; }

bool EnergyExtractor::EvaluateFrame(ctp::Topology *top) {

  std::string xmlfile = Identify() + ".xml";

  tools::Property state("state", "", "");
  tools::Property &segs = state.add("segments", "");
  tools::Property &pairs = state.add("pairs", "");
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
    if (seg->hasChrgState(+1)) {
      tools::Property &channel = segprop.add("channel", "");
      channel.setAttribute("type", "hole");
      channel.add("energy_hn",
                  (format("%1$+1.7f") % seg->getSiteEnergy(+1)).str());
      channel.add("lambda_hn",
                  (format("%1$+1.7f") % seg->getU_nC_nN(+1)).str());
      channel.add("lambda_nh",
                  (format("%1$+1.7f") % seg->getU_cN_cC(+1)).str());
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
    }
  }

  if (tools::globals::verbose) {
    // PAIRS
    ctp::QMNBList::iterator pit;
    ctp::QMNBList &nb = top->NBList();
    next = &pairs;
    for (pit = nb.begin(); pit != nb.end(); ++pit) {
      ctp::QMPair *qmp = *pit;
      tools::Property &pairprop = next->add("pair", "");
      pairprop.add("id1", (format("%1$d") % qmp->Seg1()->getId()).str());
      pairprop.add("name1", (format("%1$s") % qmp->Seg1()->getName()).str());
      pairprop.add("id2", (format("%1$d") % qmp->Seg2()->getId()).str());
      pairprop.add("name2", (format("%1$s") % qmp->Seg2()->getName()).str());

      if (qmp->isPathCarrier(+1)) {
        tools::Property &channel = pairprop.add("channel", "");
        channel.setAttribute("type", "hole");
        channel.add("deltaE_h12", (format("%1$1.7f") % qmp->getdE12(+1)).str());
        channel.add("lambdaI_h12",
                    (format("%1$1.7f") % qmp->getReorg12(+1)).str());
        channel.add("lambdaI_h21",
                    (format("%1$1.7f") % qmp->getReorg21(+1)).str());
        channel.add("lambdaO_h",
                    (format("%1$1.7f") % qmp->getLambdaO(+1)).str());
      }
      if (qmp->isPathCarrier(-1)) {
        tools::Property &channel = pairprop.add("channel", "");
        channel.setAttribute("type", "electron");
        channel.add("deltaE_e12", (format("%1$1.7f") % qmp->getdE12(-1)).str());
        channel.add("lambdaI_e12",
                    (format("%1$1.7f") % qmp->getReorg12(-1)).str());
        channel.add("lambdaI_e21",
                    (format("%1$1.7f") % qmp->getReorg21(-1)).str());
        channel.add("lambdaO_e",
                    (format("%1$1.7f") % qmp->getLambdaO(-1)).str());
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

#endif  // VOTCA_XTP_ENERGYEXTRACTOR_H
