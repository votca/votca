/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#include "votca/xtp/nodedevice.h"


using namespace std;

namespace votca {
    namespace xtp {
void NodeDevice::Compute_Self_Image_Coulomb_Potential(double startz, double device_length, Eventinfo* eventinfo) {

    double coulpot = 0.0;
    double L = device_length;
      
    double sign;
    double distz_1;
    double distz_2;
      
    for (int i=0;i<eventinfo->number_short_range_images; i++) {
        if (div(i,2).rem==0) { // even generation
            sign = -1.0;
            distz_1 = i*L + 2*startz;
            distz_2 = (i+2)*L - 2*startz; 
        }
        else { // odd generation
            sign = 1.0;
            distz_1 = (i+1)*L;
            distz_2 = (i+1)*L;
        }
        coulpot += sign*1.0/distz_1;
        coulpot += sign*1.0/distz_2;
    }
    _self_image = coulpot;    
}


        
    }
}
