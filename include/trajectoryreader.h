/// \addtogroup csg
///@{
// 
// File:   trajectoryreader.h
// Author: ruehle
//
// Created on April 20, 2007, 12:30 PM
//

#ifndef _trajectoryreader_H
#define	_trajectoryreader_H

#include <string>
#include "configuration.h"
#include "fileformatfactory.h"

using namespace std;

/**
    \brief trajectoryreader interface
    
    This class defines the interface a trajectory reader has to implement
 */
class TrajectoryReader
{
    public:
        virtual ~TrajectoryReader() {}
        /// open a trejectory file
        virtual bool Open(const string &file) = 0;
        
        virtual void Close() {};
        
        /// read in the first frame
        virtual bool FirstFrame(Configuration &conf) = 0;
        /// read in the next frame
        virtual bool NextFrame(Configuration &conf) = 0;

        virtual void RegisteredAt(ObjectFactory<string, TrajectoryReader> &factory) {}    

        static void RegisterPlugins(void);
        
        virtual TrajectoryReader *Clone(void) = 0;
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryReader> &TrjReaderFactory()
{
    static FileFormatFactory<TrajectoryReader> _TrjReaderFactory;
    return _TrjReaderFactory;
}

#endif	/* _trajectoryreader_H */

/// @}
