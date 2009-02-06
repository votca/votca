/* 
 * File:   trajectorywriter.h
 * Author: ruehle
 *
 * Created on September 7, 2007, 3:01 PM
 */

#ifndef _trajectorywriter_H
#define	_trajectorywriter_H

#include <string>
#include "fileformatfactory.h"
#include "topology.h"

using namespace std;

class TrajectoryWriter
{
public:
    TrajectoryWriter() {}
    virtual ~TrajectoryWriter() {}
    
    virtual void Open(string file, bool bAppend = false) {}
    virtual void Close() {};
    
    virtual void Write(Topology *top) {}

    virtual void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory) {}    
    
    static void RegisterPlugins(void);

    virtual TrajectoryWriter *Clone(void) = 0;
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryWriter> &TrjWriterFactory()
{
    static FileFormatFactory<TrajectoryWriter> _TrjWriterFactory;
    return _TrjWriterFactory;
}

#endif	/* _trajectorywriter_H */

