/* 
 * File:   beadtype.h
 * Author: ruehle
 *
 * Created on November 26, 2008, 12:44 PM
 */

#ifndef _BEADTYPE_H
#define	_BEADTYPE_H

#include <string>
#include "topologyitem.h"

using namespace std;

class  BeadType : public TopologyItem {
public:    
    const int &getId() const { return _id; }
    const string &getName() const { return _name; }
    
private:
    int _id;
    string _name;
    
    BeadType(Topology *parent, int id, const string &name)
    : _id(id), _name(name), TopologyItem(parent) {}
    friend class Topology;
};

#endif	/* _BEADTYPE_H */

