// 
// File:   nematicorder.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:09 PM
//

#include "nematicorder.h"
#include <tools/matrix.h>


void NematicOrder::Process(Configuration &conf)
{
    _mu.ZeroMatrix();
    _mv.ZeroMatrix();
    _mw.ZeroMatrix();
    int N=0;
    for(int n=0; n<conf.getTopology()->BeadCount();++n) {
        if( conf.getTopology()->getBead(n)->getSymmetry() ==1 )
            continue;
            
        if(conf.HasU()) {
            _mu += conf.getU(n)|conf.getU(n);
            _mu[0][0] -= 1./3.;
            _mu[1][1] -= 1./3.;
            _mu[2][2] -= 1./3.;
        }
        
        if(conf.HasV()) {
            _mv += conf.getV(n)|conf.getV(n);
            _mv[0][0] -= 1./3.;
            _mv[1][1] -= 1./3.;
            _mv[2][2] -= 1./3.;
        }
        
        if(conf.HasW()) {
            _mw += conf.getW(n)|conf.getW(n);
            _mw[0][0] -= 1./3.;
            _mw[1][1] -= 1./3.;
            _mw[2][2] -= 1./3.;
        }
        N++;
    }
    
    double f = 1./(double)N*3./2.;
    _mu = f*_mu;_mv = f*_mv;_mw = f*_mw;
    
    if(conf.HasU())
        _mu.SolveEigensystem(_nemat_u);
    if(conf.HasV())
        _mv.SolveEigensystem(_nemat_v);
    if(conf.HasW())
        _mw.SolveEigensystem(_nemat_w);               
}
