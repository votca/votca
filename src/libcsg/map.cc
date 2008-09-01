/// \addtogroup csg
///@{
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
        
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if(in.getConfiguration()->HasPos())
            cg += (*iter)._weight * in.getBeadPos((*iter)._in);
        if(in.getConfiguration()->HasVel())
            vel += (*iter)._weight * in.getBeadVel((*iter)._in);
        if(in.getConfiguration()->HasF())
            f += (*iter)._weight * in.getBeadF((*iter)._in);
    }
    if(out.getConfiguration()->HasPos())
        out.BeadPos(_out) = cg;
    if(out.getConfiguration()->HasVel())
        out.BeadVel(_out) = vel;
    if(out.getConfiguration()->HasF())
        out.BeadF(_out) = f;
}

/// \todo implement this function
void Map_Ellipsoid::Apply(Molecule &in, Molecule &out)
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), c(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    matrix m(0.);
        
    int n;
    n = 0;
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if(in.getConfiguration()->HasPos())
            cg += (*iter)._weight * in.getBeadPos((*iter)._in);
        if(in.getConfiguration()->HasVel())
            vel += (*iter)._weight * in.getBeadVel((*iter)._in);
        if(in.getConfiguration()->HasF())
            f += (*iter)._weight * in.getBeadF((*iter)._in);
         
        if((*iter)._weight>0) {
            c += in.getBeadPos((*iter)._in);
            n++;
        }
    }
    
    if(out.getConfiguration()->HasPos())
        out.BeadPos(_out) = cg;
    if(out.getConfiguration()->HasVel())
        out.BeadVel(_out) = vel;
    if(out.getConfiguration()->HasF())
        out.BeadF(_out) = f;
    
    // calculate the tensor of gyration
    c=c/(double)n;    
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if((*iter)._weight == 0) continue;
        vec v = (in.getBeadPos((*iter)._in) - c);
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
    vec v = in.getBeadPos(_matrix[1]._in) - in.getBeadPos(_matrix[0]._in);
    v.normalize();
    
    out.getConfiguration()->HasV(true);
    out.BeadV(_out) = v;    
    
    vec w = in.getBeadPos(_matrix[2]._in) - in.getBeadPos(_matrix[0]._in);
    w.normalize();
    if((v^w)*u < 0) u=vec(0.,0.,0.)-u;
    out.getConfiguration()->HasU(true);
    out.BeadU(_out) = u;
    //out.BeadV(_out) = v;
    
    //out.BeadW(_out) = es.eigenvecs[2];
}
/// @}
