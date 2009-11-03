// 
// File:   gmxtrajectorywriter.cc
// Author: ruehle
//
// Created on August 31, 2007, 2:25 PM
//


#include <string>
#include "gmxtrajectorywriter.h"

void GMXTrajectoryWriter::Open(string file, bool bAppend)
{
    //char c[1] = bAppend ? "a" : "w";
    _file = gmx::open_trx((char *)file.c_str(), "w");
}

void GMXTrajectoryWriter::Close()
{
    gmx::close_trx(_file);
}

void GMXTrajectoryWriter::Write(Topology *conf)
{
    static int step=0;   
    int N = conf->BeadCount();
    gmx::t_trxframe frame;
    gmx::rvec *x = new gmx::rvec[N];
    matrix box = conf->getBox();
    
    frame.natoms = N;
    frame.bTime = true;
    frame.time = conf->getTime();
    frame.bStep = true;
    frame.step = conf->getStep();;
    frame.x = x;
    frame.bTitle=false;
    frame.bLambda=false;
    frame.bAtoms=false;
    frame.bPrec=false;
    frame.bV=false;
    frame.bF=false;
    frame.bBox=true;
    
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            frame.box[i][j] = box[i][j];
    
    
for(int i=0; i<N; ++i) {
        vec v = conf->getBead(i)->getPos();
        x[i][0] = v.getX();
        x[i][1] = v.getY();
        x[i][2] = v.getZ(); 
    }
        
#ifdef GMX4CVS
    gmx::write_trxframe(_file, &frame, NULL);
#else
    gmx::write_trxframe(_file, &frame);
#endif
    
    step++;
    delete[] x;
}
