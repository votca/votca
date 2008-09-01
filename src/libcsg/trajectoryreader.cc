/// \addtogroup csg
///@{
// 
// File:   trajectoryreader.cc
// Author: victor
//
// Created on 27. Dezember 2007, 20:27
//


#include "trajectoryreader.h"
#include "modules/io/gmxtrajectoryreader.h"
#include "modules/io/mdptrajectoryreader.h"

void TrajectoryReader::RegisterPlugins(void)
{
    TrjReaderFactory().Register("trr", new GMXTrajectoryReader(), false);
    TrjReaderFactory().Register("xtc", new GMXTrajectoryReader(), false);
    TrjReaderFactory().Register("gro", new GMXTrajectoryReader(), false);  // not tested!
    TrjReaderFactory().Register("pdb", new GMXTrajectoryReader(), false);  // not tested!
    TrjReaderFactory().Register("mdp", new MDPTrajectoryReader(), false);
}
/// @}
