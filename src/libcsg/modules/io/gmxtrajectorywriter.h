// 
// File:   gmxtrajectorywriter.h
// Author: victor
//
// Created on 27. Dezember 2007, 20:44
//

#ifndef _GMXTRAJECTORYWRITER_H
#define	_GMXTRAJECTORYWRITER_H

#include "configuration.h"
#include "topology.h"
#include "trajectorywriter.h"

namespace gmx {
   extern "C" {
        #include <gromacs/statutil.h>
        #include <gromacs/typedefs.h>
        #include <gromacs/smalloc.h>
        #include <gromacs/vec.h>
        #include <gromacs/copyrite.h>
        #include <gromacs/statutil.h>
        #include <gromacs/tpxio.h>
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
    void Write(Configuration *conf);

    void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory);

    TrajectoryWriter *Clone() { return new GMXTrajectoryWriter(); }

    private:
        int _file;
};

#endif	/* _GMXTRAJECTORYWRITER_H */

