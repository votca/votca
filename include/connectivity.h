// 
// File:   connectivity.h
// Author: ruehle
//
// Created on January 18, 2008, 5:06 PM
//

#ifndef _CONNECTIVITY_H
#define	_CONNECTIVITY_H

#include <votca/tools/vec.h>
#include "interaction.h"

/**
    \brief Calculates connectivity of topology

    This class calculates the connectivity of the topology. It parses
    all the bonded interactions and creates a connection list.
 
 */
class Connectivity
{
public:
    struct partner_t {
        int _bead;     // the id of the neighbour
        Interaction *_inter;
    };
    
    typedef list<partner_t> container; 
    struct entry_t {
        container _neighbours;
    };
protected:
    vector<entry_t *> _list;
}


#endif	/* _CONNECTIVITY_H */

