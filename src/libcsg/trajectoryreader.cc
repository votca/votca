// 
// File:   trajectoryreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 20:27
//


#include "trajectoryreader.h"
#include "modules/io/gmxtrajectoryreader.h"

void TrajectoryReader::RegisterPlugins(void)
{
    TrjReaderFactory().Register<GMXTrajectoryReader>("trr");
    TrjReaderFactory().Register<GMXTrajectoryReader>("xtc");
    TrjReaderFactory().Register<GMXTrajectoryReader>("gro");
    TrjReaderFactory().Register<GMXTrajectoryReader>("pdb");
}
