/// \addtogroup csg
///@{
// 
// File:   trajectorywriter.cc
// Author: victor
//
// Created on 27. Dezember 2007, 20:28
//

#include <iostream>
#include "trajectorywriter.h"
#include "modules/io/pdbwriter.h"
#include "modules/io/gmxtrajectorywriter.h"
#include "modules/io/growriter.h"

using namespace std;

void TrajectoryWriter::RegisterPlugins()
{
    TrjWriterFactory().Register("pdb", new PDBWriter(), false);
    TrjWriterFactory().Register("trr", new GMXTrajectoryWriter(), false);
    TrjWriterFactory().Register("xtc", new GMXTrajectoryWriter(), false);
    TrjWriterFactory().Register("gro", new GROWriter(), false);
}
/// @}
