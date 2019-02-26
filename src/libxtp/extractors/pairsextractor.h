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

#ifndef VOTCA_XTP_PAIRSEXTRACTOR_H
#define VOTCA_XTP_PAIRSEXTRACTOR_H

#include <boost/format.hpp>
#include <votca/ctp/qmcalculator.h>
#include <votca/tools/propertyiomanipulator.h>

namespace votca {
namespace xtp {

class PairsExtractor : public ctp::QMCalculator {
 public:
  PairsExtractor(){};
  ~PairsExtractor(){};

  std::string Identify() { return "extract.pairs"; }
  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
};

void PairsExtractor::Initialize(tools::Property *options) { return; }

bool PairsExtractor::EvaluateFrame(ctp::Topology *top) {

  // Rigidify std::system (if possible)
  if (!top->Rigidify()) return 0;

  std::string xmlfile = Identify() + ".xml";

  tools::Property state("state", "", "");
  tools::Property &pairs = state.add("pairs", "");

  using boost::format;

  // PAIRS
  ctp::QMNBList::iterator pit;
  ctp::QMNBList &nb = top->NBList();
  for (pit = nb.begin(); pit != nb.end(); ++pit) {
    ctp::QMPair *qmp = *pit;
    tools::Property &pairprop = pairs.add("pair", "");
    pairprop.add("id1", (format("%1$d") % qmp->Seg1()->getId()).str());
    pairprop.add("name1", (format("%1$s") % qmp->Seg1()->getName()).str());
    pairprop.add("id2", (format("%1$d") % qmp->Seg2()->getId()).str());
    pairprop.add("name2", (format("%1$s") % qmp->Seg2()->getName()).str());

    pairprop.add("pairvector",
                 (format("%1$+1.4f %2$+1.4f %3$+1.4f") % qmp->R().getX() %
                  qmp->R().getY() % qmp->R().getZ())
                     .str());

    if (qmp->isPathCarrier(+1)) {
      tools::Property &channel = pairprop.add("channel", "");
      channel.setAttribute("type", "hole");
      channel.add("jeff2_h", (format("%1$1.7e") % qmp->getJeff2(+1)).str());
      channel.add("deltaE_h12", (format("%1$1.7f") % qmp->getdE12(+1)).str());
      channel.add("lambdaI_h12",
                  (format("%1$1.7f") % qmp->getReorg12(+1)).str());
      channel.add("lambdaI_h21",
                  (format("%1$1.7f") % qmp->getReorg21(+1)).str());
      channel.add("lambdaO_h", (format("%1$1.7f") % qmp->getLambdaO(+1)).str());
      channel.add("rate_h12", (format("%1$1.7e") % qmp->getRate12(+1)).str());
      channel.add("rate_h21", (format("%1$1.7e") % qmp->getRate12(+1)).str());
    }
    if (qmp->isPathCarrier(-1)) {
      tools::Property &channel = pairprop.add("channel", "");
      channel.setAttribute("type", "electron");
      channel.add("jeff2_e", (format("%1$1.7e") % qmp->getJeff2(-1)).str());
      channel.add("deltaE_e12", (format("%1$1.7f") % qmp->getdE12(-1)).str());
      channel.add("lambdaI_e12",
                  (format("%1$1.7f") % qmp->getReorg12(-1)).str());
      channel.add("lambdaI_e21",
                  (format("%1$1.7f") % qmp->getReorg21(-1)).str());
      channel.add("lambdaO_e", (format("%1$1.7f") % qmp->getLambdaO(-1)).str());
      channel.add("rate_e12", (format("%1$1.7e") % qmp->getRate12(-1)).str());
      channel.add("rate_e21", (format("%1$1.7e") % qmp->getRate12(-1)).str());
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

#endif  // VOTCA_XTP_PAIRSEXTRACTOR_H
