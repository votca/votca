/// \addtogroup csg
///@{
// 
// File:   gmxtrajectoryreader.cc
// Author: ruehle
//
// Created on April 5, 2007, 2:42 PM
//

#include <iostream>
#include "topology.h"
#include "gmxtrajectoryreader.h"

using namespace std;

bool GMXTrajectoryReader::Open(const string &file)
{
       _filename = file;
       return true;
}

void GMXTrajectoryReader::Close()
{
    gmx::close_trx(_gmx_status);
}

bool GMXTrajectoryReader::FirstFrame(Configuration &conf)
{
    if(!gmx::read_first_frame(&_gmx_status,(char*)_filename.c_str(),&_gmx_frame,TRX_READ_X | TRX_READ_F))
        return false;
    matrix m;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            m[i][j] = _gmx_frame.box[j][i];
    conf.setBox(m);
    conf.setTime(_gmx_frame.time);
    conf.setStep(_gmx_frame.step);
    cout << endl;
    
    conf.HasPos(true);
    conf.HasF(_gmx_frame.bF);
    
    for(int i=0; i<_gmx_frame.natoms; i++) {
        double r[3] = { _gmx_frame.x[i][XX],  _gmx_frame.x[i][YY], _gmx_frame.x[i][ZZ] };
        //CBeadInfo *bi = conf.getTopology()->getBeads()[i];
        conf.setPos(i, r);
        if(_gmx_frame.bF) {
            double f[3] = { _gmx_frame.f[i][XX],  _gmx_frame.f[i][YY], _gmx_frame.f[i][ZZ] };        
            conf.setF(i, f);        
        }
    }
    return true;
}

bool GMXTrajectoryReader::NextFrame(Configuration &conf)
{
    if(!gmx::read_next_frame(_gmx_status,&_gmx_frame))
        return false;
    matrix m;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            m[i][j] = _gmx_frame.box[j][i];
    conf.setTime(_gmx_frame.time);
    conf.setStep(_gmx_frame.step);
    conf.setBox(m);
    
    conf.HasF(_gmx_frame.bF);
    
    for(int i=0; i<_gmx_frame.natoms; i++) {
        double r[3] = { _gmx_frame.x[i][XX],  _gmx_frame.x[i][YY], _gmx_frame.x[i][ZZ] };
        conf.setPos(i, r);
        if(_gmx_frame.bF) {
            double f[3] = { _gmx_frame.f[i][XX],  _gmx_frame.f[i][YY], _gmx_frame.f[i][ZZ] };        
            conf.setF(i, f);        
        }
    }
    return true;
}
/// @}
