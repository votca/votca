/* 
 * File:   topologyitem.h
 * Author: ruehle
 *
 * Created on January 27, 2009, 5:02 PM
 */

#ifndef _TOPOLOGYITEM_H
#define	_TOPOLOGYITEM_H

class Topology;

class TopologyItem
{
protected:
    TopologyItem(Topology *parent)
        : _parent(parent) {}
    
    Topology *_parent;
    
    friend class Topology;
};

#endif	/* _TOPOLOGYITEM_H */

