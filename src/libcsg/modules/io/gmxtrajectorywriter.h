// 
// File:   gmxtrajectorywriter.h
// Author: victor
//
// Created on 27. Dezember 2007, 20:44
//

#ifndef _GMXTRAJECTORYWRITER_H
#define	_GMXTRAJECTORYWRITER_H

#include "topology.h"
#include "trajectorywriter.h"

namespace gmx {
   extern "C" {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool
}
using namespace std;

class GMXTrajectoryWriter
    : public TrajectoryWriter
{
public:
    
    void Open(string file, bool bAppend = false);
    void Close();
    void Write(Topology *conf);

    void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory);

    TrajectoryWriter *Clone() { return new GMXTrajectoryWriter(); }

    private:
        int _file;
};

#endif	/* _GMXTRAJECTORYWRITER_H */

