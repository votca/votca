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

namespace votca { namespace kmc {
  
using namespace std;

enum EventType { Injection, Collection, Transfer, Recombination, Dissociation, ExcitonFormation };


class Event {
    
public:
    
    Event(EventType type, Carrier *carrier1, Carrier *carrier2) : _type(type), _carrier1(carrier1), _carrier2(carrier2) { ; }
    
    virtual void OnExecute() { ; }
    const double &getRate() const { return _rate; }
    
    
private:
    
    EventType _type;
    double _rate;
    Carrier *_carrier1;
    Carrier *_carrier2;   

    
};

class Recombination : Event {
public:
    
    void OnExecute() { ; }
    // TODO updating the state
    
private:    
};


}} 

#endif

