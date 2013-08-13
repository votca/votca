/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_EVENT_H_
#define __VOTCA_KMC_EVENT_H_

#include <votca/kmc/carrier.h>

namespace votca { namespace kmc {
  
using namespace std;

enum EventType { _ElectronTransfer, _HoleTransfer, _ElectronInjection, _HoleInjection, _ElectronCollection, _HoleCollection, _Recombination, _Dissociation};

/*
 * Abstract base class for all events  
 */
class Event {
    
public:
    
    void Initialize(Carrier *carrier1, Carrier *carrier2) { ; }
    
    virtual void onExecute() = 0;
    
    virtual EventType &getEventType(){ return _type; }
    
    // maybe make _rate public (check the efficiency)) or inline it
    const double &getRate() const { return _rate; }
    
    
protected:
    
    EventType _type;
    double _rate;
    Carrier *_carrier1; 
    Carrier *_carrier2;
    
};



}} 

#endif

