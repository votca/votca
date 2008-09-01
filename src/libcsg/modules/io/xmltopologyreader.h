/// \addtogroup csg
///@{
// 
// File:   pdbtopologyreader.h
// Author: victor
//
// Created on 14. April 2008, 13:33
//

#ifndef _XMLTOPOLOGYREADER_H
#define	_XMLTOPOLOGYREADER_H

#include <string>
#include <libxml/xmlreader.h>
#include "topologyreader.h"
#include <stack>

using namespace std;
    
/**

*/
class XMLTopologyReader
   : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(string file, Topology &top);

    TopologyReader *Clone() { return new XMLTopologyReader(); }

private:    
    void ReadTopolFile(string file);
          
    void ParseTopology(xmlNodePtr node);
    void ParseMolecules(xmlNodePtr node);
    
    Topology *_top;
};

#endif	/* _PDBTOPOLOGYREADER_H */

/// @}
