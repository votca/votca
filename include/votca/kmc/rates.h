/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_RATES_H_
#define __VOTCA_KMC_RATES_H_

#include <votca/kmc/bsumtree.h>

namespace votca { namespace kmc {
  
using namespace std;

class Rates { // Bookkeeping of rates
public:
    
    void initialize(int injection_events, int ncharges, int nbtot);
    void set_injection_rate(int injection_event, double rate);
    void set_jump_rate(int charge, int jump, double rate);
    double recompute();
    void search(int& charge, int& jump, double searchkey);
    void resize(int ncharges);
    
private:
  Bsumtree sumtree;
  int injection_events;
  long ncharges;
  int nbtot;
};

void Rates::initialize(int injection_events_, int ncharges_, int nbtot_) { // Must be called before use
  injection_events = injection_events_;
  ncharges = ncharges_;
  nbtot = nbtot_;
  sumtree.initialize(injection_events + ncharges*nbtot);
}

void Rates::set_injection_rate(int injection_event, double rate) {
  sumtree.setelement(injection_event,rate);
}
  
void Rates::set_jump_rate(int charge, int jump, double rate) {
  int i = injection_events + charge*nbtot+jump;
  sumtree.setelement(i,rate);
}

double Rates::recompute() {
  return sumtree.compute_sum(); // Calculate total rate
}

void Rates::search(int& charge, int& jump, double searchkey) {
  long i = sumtree.search(searchkey);
  if (i<injection_events) { // Injection occurred
    charge = -1; // Encoding for an injection event
    jump = i; // Encoding for which injection event occurred
  } 
  else {
    i -= injection_events;
    // Reconstruct charge and jump numbers
    ldiv_t mydiv = ldiv(i,(long)nbtot);
    charge = mydiv.quot;
    jump = mydiv.rem;
  }
}

void Rates::resize(int ncharges_) { // Resize sumtree to accomodate more charges
  ncharges = ncharges_;
  sumtree.resize(injection_events + ncharges*nbtot);
}

}}

#endif