// 
// File:   nematicorder.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:09 PM
//

#include "nematicorder.h"
#include <tools/matrix.h>
#include "topology.h"

void NematicOrder::Process(Topology &top)
{
    _mu.ZeroMatrix();
    _mv.ZeroMatrix();
    _mw.ZeroMatrix();
    int N=0;
    bool bU, bV, bW;
    bU=bV=bW = false;
    
    for(BeadContainer::iterator iter = top.Beads().begin();
    iter!=top.Beads().end();++iter) {
        Bead *bead = *iter;
        
        if( bead->getSymmetry() ==1 )
            continue;
            
        if(bead->HasU()) {
            _mu += bead->getU()|bead->getU();
            _mu[0][0] -= 1./3.;
            _mu[1][1] -= 1./3.;
            _mu[2][2] -= 1./3.;
            bU = true;
        }
        
        if(bead->HasV()) {
            _mv += bead->getV()|bead->getV();
            _mv[0][0] -= 1./3.;
            _mv[1][1] -= 1./3.;
            _mv[2][2] -= 1./3.;
            bV = true;
        }
        
        if(bead->HasW()) {
            _mw += bead->getW()|bead->getW();
            _mw[0][0] -= 1./3.;
            _mw[1][1] -= 1./3.;
            _mw[2][2] -= 1./3.;
            bU = false;
        }
        N++;
    }
    
    double f = 1./(double)N*3./2.;
    _mu = f*_mu;_mv = f*_mv;_mw = f*_mw;
    
    if(bU)
        _mu.SolveEigensystem(_nemat_u);
    if(bV)
        _mv.SolveEigensystem(_nemat_v);
    if(bW)
        _mw.SolveEigensystem(_nemat_w);               
}
