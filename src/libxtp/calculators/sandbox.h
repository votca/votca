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


#ifndef SANDBOX2_H
#define SANDBOX2_H

#include <votca/xtp/parallelpaircalc.h>

namespace votca { namespace xtp {

class Sandbox : public ParallelPairCalculator
{

public:
    
    Sandbox() { };
   ~Sandbox() { };

    string  Identify() { return "Sandbox"; }
    void    Initialize(Property *options);
    using ParallelPairCalculator::EvalPair;
    void    EvalPair(Topology *top, QMPair *qmpair, int slot);

};


void Sandbox::Initialize(Property *options) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );

    _nThreads = 1;

    cout << endl << "... ... Initialize with " << _nThreads << " threads ";    

}


void Sandbox::EvalPair(Topology *top, QMPair *qmpair, int slot) {
  
    this->LockCout();
    cout << "\r... ... Overloading pair " << qmpair->getId() << ". " << flush;
    this->UnlockCout();

    //int ij;
    for (int i = 0; i < 2000; i++) {
        for (int j = 0; j < 2000; j++) {
            //ij = i+j;
        }
    }

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin(); sit != top->Segments().end(); sit++) {
        //Segment *seg = *sit;

        vector<Atom*> ::iterator ait;
        for (ait= top->Atoms().begin(); ait != top->Atoms().end(); ait++) {
            //Atom *atm = *ait;

            //int id = atm->getId();

        }

        //Atom *atm = seg->Atoms()[1];
        //atm = seg->Atoms()[2];
        //atm = seg->Atoms()[10000];
        //atm = seg->Atoms()[14300];

    }


}













}} /* exit namespace votca::xtp */

#endif /* SANDBOX2_H */
