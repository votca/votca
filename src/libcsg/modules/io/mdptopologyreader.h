/// \addtogroup csg
///@{
/* 
 * File:   mdptopologyreader.h
 * Author: ruehle
 *
 * Created on June 12, 2008, 10:47 AM
 */

#ifndef _MDPTOPOLOGYREADER_H
#define	_MDPTOPOLOGYREADER_H

#include "topologyreader.h"

/**
    \brief reader for Alexander Lyubartsev's md format
 
 */
class MDPTopologyReader
    : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);
    TopologyReader *Clone() { return new MDPTopologyReader(); }
    
private:
       
};


#endif	/* _MDPTOPOLOGYREADER_H */

/// @}
