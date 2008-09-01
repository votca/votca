/// \addtogroup csg
///@{
// 
// File:   gmxtopologyreader.h
// Author: ruehle
//
// Created on April 5, 2007, 12:33 PM
//

#ifndef _gmxtopologyreader_H
#define	_gmxtopologyreader_H

#include <string>
#include "topologyreader.h"

using namespace std;
    
/**
    \brief reader for gromacs topology files

    This class encapsulates the gromacs reading functions and provides an interface to fill a topolgy class

*/
class GMXTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);
    
    TopologyReader *Clone() { return new GMXTopologyReader(); }

private:
};

#endif	/* _gmxtopologyreader_H */

/// @}
