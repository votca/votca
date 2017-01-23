/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef _PAIRCALCULATOR2_H
#define	_PAIRCALCULATOR2_H

#include <votca/xtp/qmcalculator.h>

namespace votca { namespace ctp {

class PairCalculator2 : public XQMCalculator
{
public:

    PairCalculator2() {};
    virtual ~PairCalculator2() {};

    bool EvaluateFrame(Topology *top);
    virtual void EvaluatePair(Topology *top, QMPair *pair) { };
};

bool PairCalculator2::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;

    QMNBList &nblist = top->NBList();

    QMNBList::iterator pit;
    for (pit = nblist.begin(); pit != nblist.end(); pit++) {

        EvaluatePair(top, *pit);

        if ( (*pit)->getId() == -1 ) {
            
            string pairId = boost::lexical_cast<string>((*pit)->getId());
            string pdbname = "Pair" + pairId + ".pdb";
            (*pit)->WritePDB(pdbname);
        }

    }
    
    return 1;
}

}}

#endif	/* _PAIRCALCULATOR2_H */

 
