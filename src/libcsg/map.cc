// 
// File:   map.cc
// Author: ruehle
//
// Created on April 13, 2007, 4:41 PM
//

#include "map.h"
#include <iostream>
#include <tools/matrix.h>

using namespace std;

Map::~Map()
{
    vector<BeadMap *>::iterator iter;
    for(iter=_maps.begin();iter!=_maps.end();++iter)
        delete (*iter);    
    _maps.clear();
}

void Map::Apply(Molecule &in, Molecule &out)
{
    vector<BeadMap *>::iterator iter;
    for(iter=_maps.begin();iter!=_maps.end();++iter)
        (*iter)->Apply(in, out);
}

void Map_Sphere::Apply(Molecule &in, Molecule &out)
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    bool bPos, bVel, bF;
    bPos=bVel=bF=false;
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if(in.getBead((*iter)._in)->HasPos()) {
            cg += (*iter)._weight * in.getBead((*iter)._in)->getPos();
            bPos=true;
        }
        if(in.getBead((*iter)._in)->HasVel()) {
            vel += (*iter)._weight * in.getBead((*iter)._in)->getVel();
            bVel = true;
        }
        if(in.getBead((*iter)._in)->HasF()) {
            /// \todo fix me, right calculation should be F_i = m_cg / sum(w_i) * sum(w_i/m_i*F_i)
            //f += (*iter)._weight * in.getBeadF((*iter)._in);
            f += in.getBead((*iter)._in)->getF();
            bF = true;
        }
    }
    if(bPos)
        out.getBead(_out)->setPos(cg);
    if(bVel)
        out.getBead(_out)->setVel(vel);
    if(bF)
        out.getBead(_out)->setF(f);
}

/// \todo implement this function
void Map_Ellipsoid::Apply(Molecule &in, Molecule &out)
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), c(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    matrix m(0.);
     bool bPos, bVel, bF;
    bPos=bVel=bF=false;
      
    int n;
    n = 0;
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if(in.getBead((*iter)._in)->HasPos()) {
            cg += (*iter)._weight * in.getBead((*iter)._in)->getPos();
            bPos=true;
        }
        if(in.getBead((*iter)._in)->HasVel() == true) {
            vel += (*iter)._weight * in.getBead((*iter)._in)->getVel();
            bVel = true;
        }
        if(in.getBead((*iter)._in)->HasF()) {
            /// \todo fix me, right calculation should be F_i = m_cg / sum(w_i) * sum(w_i/m_i*F_i)
            //f += (*iter)._weight * in.getBeadF((*iter)._in);
            f += in.getBead((*iter)._in)->getF();
            bF = true;
        }
        
        if((*iter)._weight>0) {
            c += in.getBead((*iter)._in)->getPos();
            n++;
        }
    }
    
    if(bPos)
        out.getBead(_out)->setPos(cg);
    if(bVel)
        out.getBead(_out)->setVel(vel);
    if(bF)
        out.getBead(_out)->setF(f);
    
    // calculate the tensor of gyration
    c=c/(double)n;    
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if((*iter)._weight == 0) continue;
        vec v = (in.getBead((*iter)._in)->getPos() - c);
        //v = vec(1, 0.5, 0) * 0.*(drand48()-0.5)
        //    + vec(0.5, -1, 0) * (drand48()-0.5)
        //    + vec(0, 0, 1) * (drand48()-0.5);
        m[0][0] += v.getX()*v.getX();
        m[0][1] += v.getX()*v.getY();
        m[0][2] += v.getX()*v.getZ();        
        m[1][1] += v.getY()*v.getY();
        m[1][2] += v.getY()*v.getZ();
        m[2][2] += v.getZ()*v.getZ();
    }
    m[1][0] = m[0][1];
    m[2][0] = m[0][2];
    m[2][1] = m[1][2];
    
    // calculate the eigenvectors
    matrix::eigensystem_t es;
    m.SolveEigensystem(es);
    
    vec u = es.eigenvecs[0];
    vec v = in.getBead(_matrix[1]._in)->getPos() - in.getBead(_matrix[0]._in)->getPos();
    v.normalize();
    
    out.getBead(_out)->setV(v);    
    
    vec w = in.getBead(_matrix[2]._in)->getPos() - in.getBead(_matrix[0]._in)->getPos();
    w.normalize();
    
    
    
    if((v^w)*u < 0) u=vec(0.,0.,0.)-u;
    out.getBead(_out)->setU(u);
    
    //write out w
    w=u^v;
    w.normalize();
    out.getBead(_out)->setW(w);
    
    //out.BeadV(_out) = v;
    
    //out.BeadW(_out) = es.eigenvecs[2];
}
