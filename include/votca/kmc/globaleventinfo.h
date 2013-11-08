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

#ifndef __VOTCA_KMC_GLOBALEVENTINFO_H_
#define __VOTCA_KMC_GLOBALEVENTINFO_H_

namespace votca { namespace kmc {
  
using namespace std;

class Globaleventinfo {

public:
    
    double alpha;
    double beta;
    double efield;
    double injection_barrier;
    double binding_energy;
    double coulomb_strength;
    double coulcut;
    double self_image_prefactor;
    
    bool left_injection[2];
    bool right_injection[2];
    bool device;
    string formalism;
    
    int nr_sr_images;
    long nr_of_lr_images;
    int state_grow_size;
    
    double electron_prefactor;
    double hole_prefactor;
    double injection_prefactor;
    double recombination_prefactor;
    double collection_prefactor;
};

}} 

#endif

