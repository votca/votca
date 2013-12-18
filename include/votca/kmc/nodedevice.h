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

#ifndef _VOTCA_KMC_NODEDEVICE_H
#define	_VOTCA_KMC_NODEDEVICE_H

#include <votca/kmc/nodesql.h>
#include <votca/kmc/eventinfo.h>

namespace votca { namespace kmc {

class NodeDevice : public NodeSQL
{
public:
    
    NodeDevice(int id, tools::vec position) : NodeSQL(id, position){};

    void setLayer(int layer){_layer = layer;}
//    void Compute_Self_Coulomb_potential(double startx, votca::tools::vec simboxsize, Eventinfo* eventinfo); 
    void Compute_Self_Image_Coulomb_Potential(double startx, double device_length, Eventinfo* eventinfo);
    void setSelfImage(double self_image) {_self_image = self_image;}
    
    const double &self_image() const { return _self_image; }     
    const int &layer() const {return _layer;}

    
private:

    int _layer;
    double _self_image;

    ;
    
};

void NodeDevice::Compute_Self_Image_Coulomb_Potential(double startx, double device_length, Eventinfo* eventinfo) {

    double coulpot = 0.0;
    double L = device_length;
      
    int sign;
    double distx_1;
    double distx_2;
    bool outside_cut_off1 = false;
    bool outside_cut_off2 = false;
      
    while(!(outside_cut_off1&&outside_cut_off2)) {
        for (int i=0;i<eventinfo->nr_sr_images; i++) {
            if (div(i,2).rem==0) { // even generation
                sign = -1;
                distx_1 = i*L + 2*startx;
                distx_2 = (i+2)*L - 2*startx; 
            }
            else {
                sign = 1;
                distx_1 = (i+1)*L;
                distx_2 = (i+1)*L;
            }
            if (distx_1<=eventinfo->coulcut) {
                coulpot += sign*1.0/sqrt(distx_1)-1.0/(eventinfo->coulcut);
            }
            else {
                outside_cut_off1 = true;
            }
            if (distx_2<=eventinfo->coulcut) {
                coulpot += sign*1.0/sqrt(distx_2)-1.0/(eventinfo->coulcut);
            }
            else {
                outside_cut_off2 = true;
            }
        }
    }
    _self_image = eventinfo->coulomb_strength*eventinfo->self_image_prefactor*coulpot;    
}

}}

#endif	/* _VOTCA_KMC_NODEDEVICE_H */

