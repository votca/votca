/// \addtogroup csg
///@{
// 
// File:   pdbwriter.h
// Author: victor
//
// Created on 27. Dezember 2007, 20:41
//

#ifndef _PDBWRITER_H
#define	_PDBWRITER_H

#include <stdio.h>
#include "topology.h"
#include "trajectorywriter.h"

using namespace std;

class PDBWriter
: public TrajectoryWriter
{
public:
    
    void Open(string file, bool bAppend = false);
    void Close();
    
    void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory) {}    

    void Write(Configuration *conf);

    TrajectoryWriter *Clone() { return new PDBWriter(); }

    private:
        FILE *_out;
};

#endif	/* _PDBWRITER_H */

/// @}
