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
    TrjWriterFactory().Register<PDBWriter>("pdb");
    TrjWriterFactory().Register<GMXTrajectoryWriter>("trr");
    TrjWriterFactory().Register<GMXTrajectoryWriter>("xtc");
    TrjWriterFactory().Register<GROWriter>("gro");
}
