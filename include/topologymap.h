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

