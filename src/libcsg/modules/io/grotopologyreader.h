// 
// File:   grotopologyreader.h
// Author: victor
//
// Created on 27. Dezember 2007, 21:46
//

#ifndef _GROTOPOLOGYREADER_H
#define	_GROTOPOLOGYREADER_H

#include <string>
#include "topologyreader.h"

using namespace std;
    
/**
    \brief reader for gromacs topology files

    This class encapsulates the gromacs reading functions and provides an interface to fill a topolgy class

*/
class GROTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);
    TopologyReader *Clone() { return new GROTopologyReader(); }
    
private:
};


#endif	/* _GROTOPOLOGYREADER_H */

