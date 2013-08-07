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

#ifndef __VOTCA_KMC_NODE_H_
#define __VOTCA_KMC_NODE_H_

#include <votca/tools/vec.h>

namespace votca { namespace kmc {
  
using namespace std;

enum CarrierType{ Electron, Hole, Exciton };

class Node {
    public:
        const int &getID() const { return _id; }
        const string &getType() const { return _type; }
        const vec &getPosition() const { return _position; }
//        const vector<Node*> &getNeighbours(CarrierType type){ return _neighbours[type]; }
        const vector<Node*> &getNeighbours{ return _neighbours[_id]; }
        const double &getOccupation( CarrierType type ) const { return _occupation[type]; }
        const int &getCarrierID const {return _carrierid;} // -1 if empty
        const double &getStaticEnergy const {return _static_energy;}
        
//        void setID (int index) {
//            _id = index;
//        }
//        void setPosition (double ix, double iy, double iz) {
//            _position = vec(ix,iy,iz);
//        }
        
        void setNeighbours (int vertex2) {
            _neighbours.push_back(index2);
        }
        
//        void setOccupation(CarrierType type, double occup) {
//            _occupation[type] = occup;
//        }
        
//        void setCarrierID(int carID) {
//            _carrierID = carID;
//        }
        
        ;

    private:
        int _id;
        const string _type;
        vec _position;
        vector<Node*> _neighbours;
        double _occupation[3];
        int _carrierID;
        double _static_energy;
};

}} 

#endif

