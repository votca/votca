/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
/* 
 * File:   topologymap.h
 * Author: ruehle
 *
 * Created on August 13, 2007, 3:23 PM
 */

#ifndef _topologymap_H
#define	_topologymap_H

#include "map.h"
#include "topology.h"
#include <vector>

using namespace std;

class TopologyMap 
{
public:
    ~TopologyMap();
    
    TopologyMap(Topology *in, Topology *out);
    
    void AddMoleculeMap(Map *map);

    void Apply();
    
private:
    Topology *_in;
    Topology *_out;
    
    typedef vector<Map *> MapContainer;
    MapContainer _maps;
};

inline TopologyMap::TopologyMap(Topology *in, Topology *out)
    : _in(in), _out(out) {}
    
inline void TopologyMap::AddMoleculeMap(Map *map)
{
        _maps.push_back(map);
}

#endif	/* _topologymap_H */

