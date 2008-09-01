// 
// File:   topologyreader.h
// Author: ruehle
//
// Created on January 18, 2008, 3:36 PM
//

#ifndef _TOPOLOGYREADER_H
#define	_TOPOLOGYREADER_H

#include <string>
#include "topology.h"
#include "fileformatfactory.h"

using namespace std;

class TopologyReader
{
    public:
        virtual ~TopologyReader() {}
        /// open a trejectory file
        virtual bool ReadTopology(string file, Topology &top) = 0;
        
        virtual void RegisteredAt(ObjectFactory<string, TopologyReader> &factory) {}    

        static void RegisterPlugins(void);
        
        virtual TopologyReader *Clone(void) = 0;
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TopologyReader> &TopReaderFactory()
{
    static FileFormatFactory<TopologyReader> _TopReaderFactory;
    return _TopReaderFactory;
}

#endif	/* _TOPOLOGYREADER_H */

