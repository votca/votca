// 
// File:   trajectory.h
// Author: ruehle
//
// Created on April 26, 2007, 3:55 PM
//

#ifndef _trajectory_H
#define	_trajectory_H

#include "topology.h"
#include "trajectoryreader.h"
#include "configuration.h"

class Trajectory
{
    public:
        Trajectory(Topology *top, TrajectoryReader *reader);
        
        bool FirstFrame() { return _reader->FirstFrame(_conf); }
        bool NextFrame() { return _reader->NextFrame(_conf); }
        
        Configuration &getConfiguration() { return _conf; }
    private:
        Configuration _conf;
        Topology *_topology;
        TrajectoryReader *_reader;
};

#endif	/* _trajectory_H */

