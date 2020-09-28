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
#ifndef VOTCA_XTP_PAIRCALCULATOR_H
#define VOTCA_XTP_PAIRCALCULATOR_H

// Local VOTCA includes
#include "qmcalculator.h"

namespace votca {
namespace xtp {

class PairCalculator : public QMCalculator {
 public:
  PairCalculator(){};
  virtual ~PairCalculator(){};

  bool EvaluateFrame(Topology *top);
  virtual void EvaluatePair(Topology *top, QMPair *pair){};
};

bool PairCalculator::EvaluateFrame(Topology *top) {

  // Rigidify if (a) not rigid yet (b) rigidification at all possible
  if (!top->isRigid()) {
    bool isRigid = top->Rigidify();
    if (!isRigid) {
      return 0;
    }
  } else {
    std::cout << std::endl << "... ... System is already rigidified.";
  }
  std::cout << std::endl;

  QMNBList &nblist = top->NBList();

  QMNBList::iterator pit;
  for (pit = nblist.begin(); pit != nblist.end(); pit++) {

    EvaluatePair(top, *pit);

    if ((*pit)->getId() == -1) {

      std::string pairId = boost::lexical_cast<std::string>((*pit)->getId());
      std::string pdbname = "Pair" + pairId + ".pdb";
    }
  }

  return 1;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PAIRCALCULATOR_H
