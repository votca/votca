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


#include <votca/xtp/parallelpaircalc.h>

namespace votca { namespace xtp {
    namespace CTP = votca::ctp;
bool ParallelPairCalculator::EvaluateFrame(CTP::Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { std::cout << std::endl << "... ... System is already rigidified."; }
    std::cout << std::endl;        

    std::vector<PairOperator*> pairOps;
    this->InitSlotData(top);

    _nextPair = top->NBList().begin();

    for (unsigned int id = 0; id < _nThreads; id++) {
        PairOperator *newOp = new PairOperator(id, top, this);
        pairOps.push_back(newOp);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        pairOps[id]->Start();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        pairOps[id]->WaitDone();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        delete pairOps[id];
    }

    pairOps.clear();

    this->PostProcess(top);
    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

CTP::QMPair *ParallelPairCalculator::RequestNextPair(int opId, CTP::Topology *top) {

    _nextPairMutex.Lock();

    CTP::QMPair *workOnThis;

    if (_nextPair == top->NBList().end()) {
        workOnThis = NULL;
    }
    else {
        CTP::QMPair *workOnThat = *_nextPair;
        _nextPair++;
        workOnThis = workOnThat;
    }

    _nextPairMutex.Unlock();

    return workOnThis;
}

// +++++++++++++++++++++++++++++ //
// PairOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void ParallelPairCalculator::PairOperator::Run(void) {

    while (true) {

        CTP::QMPair *qmpair = _master->RequestNextPair(_id, _top);

        if (qmpair == NULL) { break; }
        else { this->_master->EvalPair(_top, qmpair, this); }
    }
}

}}
