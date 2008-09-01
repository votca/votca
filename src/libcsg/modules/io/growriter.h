/// \addtogroup csg
///@{
// 
// File:   growriter.h
// Author: victor
//
// Created on 9. January 2008, 18:11
//

#ifndef _GROWRITER_H
#define	_GROWRITER_H

#include <stdio.h>
#include "topology.h"
#include "trajectorywriter.h"

using namespace std;

class GROWriter
: public TrajectoryWriter
{
public:
    
    void Open(string file, bool bAppend = false);
    void Close();
    
    void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory) {}    

    void Write(Configuration *conf);

    TrajectoryWriter *Clone() { return new GROWriter(); }

    private:
        FILE *_out;
};

#endif	/* _GROWRITER_H */

/// @}
