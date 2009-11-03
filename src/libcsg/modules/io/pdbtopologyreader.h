// 
// File:   pdbtopologyreader.h
// Author: victor
//
// Created on 27. Dezember 2007, 21:09
//

#ifndef _PDBTOPOLOGYREADER_H
#define	_PDBTOPOLOGYREADER_H

#include <string>
#include "topology.h"
#include "topologyreader.h"

using namespace std;
    
/**
    \brief reader for gromacs topology files

    This class encapsulates the gromacs reading functions and provides an interface to fill a topolgy class

*/
class PDBTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);

private:
};

#endif	/* _PDBTOPOLOGYREADER_H */

