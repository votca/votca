// 
// File:   map.cc
// Author: ruehle
//
// Created on April 13, 2007, 4:41 PM
//

#include "map.h"
#include <iostream>
#include <tools/matrix.h>
#include <tools/tokenizer.h>
#include <numeric>

using namespace std;

Map::~Map()
{
    vector<BeadMap *>::iterator iter;
    for(iter=_maps.begin();iter!=_maps.end();++iter)
        delete (*iter);    
    _maps.clear();
}

void Map::Apply()
{
    vector<BeadMap *>::iterator iter;
    for(iter=_maps.begin();iter!=_maps.end();++iter)
        (*iter)->Apply();
}

void Map_Sphere::Initialize(Molecule *in, Bead *out, Property *opts_bead, Property *opts_map) {
    BeadMap::Initialize(in, out, opts_bead, opts_map);
    
    vector<string> beads;
    vector<double> weights;
    vector<double> fweights;

    // get the beads
    string s(_opts_bead->get("beads").value());
    Tokenizer tok_beads(s, " \n\t");
    tok_beads.ToVector(beads);

    // get vector of weights
    Tokenizer tok_weights(_opts_map->get("weights").value(), " \n\t");
    tok_weights.ConvertToVector<double>(weights);

    // check weather weights and # beads matches
    if (beads.size() != weights.size())
        throw runtime_error(string("number of subbeads in " +
                opts_bead->get("name").value()
                + "and number of weights in map "
                + opts_map->get("name").value() + " do not match"));

    // normalize the weights
    double norm = 1./ std::accumulate(weights.begin(), weights.end(), 0.);
    
    transform(weights.begin(), weights.end(), weights.begin(), bind2nd(multiplies<double>(), norm));
    // get the d vector if exists or initialize same as weights
    vector<double> d;
    if(_opts_map->exists("d")) {
        Tokenizer tok_weights(_opts_map->get("d").value(), " \n\t");
        tok_weights.ConvertToVector(d);
        // normalize d coefficients
        norm = 1./std::accumulate(d.begin(), d.end(), 0.);
        transform(d.begin(), d.end(), d.begin(), bind2nd(multiplies<double>(), norm));
    } else {
        // initialize force-weights with weights
        d.resize(weights.size());
        copy(weights.begin(), weights.end(), d.begin());
    }

    // check weather number of d coeffs is correct
    if (beads.size() != d.size()) {
        throw runtime_error(string("number of subbeads in " +
            opts_bead->get("name").value()
            + "and number of d-coefficients in map "
            + opts_map->get("name").value() + " do not match"));
    }

    fweights.resize(weights.size());
    // calculate force weights by d_i/w_i
    for(int i=0; i<weights.size(); ++i) {
        if(weights[i] == 0 && d[i]!=0) {
            throw runtime_error(
                "A d coefficient is nonzero while weights is zero in mapping "
                + opts_map->get("name").value());
        }
        if(weights[i] != 0)
            fweights[i] = d[i] / weights[i];
        else
            fweights[i] = 0;        
    }

    for (size_t i = 0; i < beads.size(); ++i) {
        int iin = in->getBeadByName(beads[i]);
        if (iin < 0)
            throw std::runtime_error(string("mapping error: molecule " + beads[i] + " does not exist"));
        AddElem(in->getBead(iin), weights[i], fweights[i]);
    }
}

void Map_Sphere::Apply()
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    bool bPos, bVel, bF;
    bPos=bVel=bF=false;
    _out->ParentBeads().clear();
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        Bead *bead = iter->_in;
       _out->ParentBeads().push_back(bead->getId());
        if(bead->HasPos()) {
            cg += (*iter)._weight * bead->getPos();
            bPos=true;
        }
        if(bead->HasVel()) {
            vel += (*iter)._weight * bead->getVel();
            bVel = true;
        }
        if(bead->HasF()) {
            f += (*iter)._force_weight * bead->getF();
            bF = true;
        }
    }
    if(bPos)
        _out->setPos(cg);
    if(bVel)
        _out->setVel(vel);
    if(bF)
        _out->setF(f);
}

/// \todo implement this function
void Map_Ellipsoid::Apply()
{
    vector<element_t>::iterator iter;
    vec cg(0., 0., 0.), c(0., 0., 0.), f(0.,0.,0.), vel(0.,0.,0.);
    matrix m(0.);
     bool bPos, bVel, bF;
    bPos=bVel=bF=false;
      
    int n;
    n = 0;
    _out->ParentBeads().clear();
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
       Bead *bead = iter->_in;
       _out->ParentBeads().push_back(bead->getId());
       if(bead->HasPos()) {
            cg += (*iter)._weight * bead->getPos();
            bPos=true;
        }
        if(bead->HasVel() == true) {
            vel += (*iter)._weight * bead->getVel();
            bVel = true;
        }
        if(bead->HasF()) {
            /// \todo fix me, right calculation should be F_i = m_cg / sum(w_i) * sum(w_i/m_i*F_i)
            //f += (*iter)._weight * _in->getBeadF((*iter)._in);
            f += (*iter)._force_weight * bead->getF();
            bF = true;
        }
        
        if((*iter)._weight>0 && bead->HasPos()) {
            c += bead->getPos();
            n++;
        }
    }
    
    if(bPos)
        _out->setPos(cg);
    if(bVel)
        _out->setVel(vel);
    if(bF)
        _out->setF(f);
    
    // calculate the tensor of gyration
    c=c/(double)n;    
    for(iter = _matrix.begin(); iter != _matrix.end(); ++iter) {
        if((*iter)._weight == 0) continue;
        Bead *bead = iter->_in;
        vec v = bead->getPos() - c;
        //v = vec(1, 0.5, 0) * 0.*(drand48()-0.5)
        //    + vec(0.5, -1, 0) * (drand48()-0.5)
        //    + vec(0, 0, 1) * (drand48()-0.5);
        
        //Normalize the tensor with 1/number_of_atoms_per_bead
        m[0][0] += v.getX()*v.getX()/(double)_matrix.size();
        m[0][1] += v.getX()*v.getY()/(double)_matrix.size();
        m[0][2] += v.getX()*v.getZ()/(double)_matrix.size();        
        m[1][1] += v.getY()*v.getY()/(double)_matrix.size();
        m[1][2] += v.getY()*v.getZ()/(double)_matrix.size();
        m[2][2] += v.getZ()*v.getZ()/(double)_matrix.size();
    }
    m[1][0] = m[0][1];
    m[2][0] = m[0][2];
    m[2][1] = m[1][2];
    
    // calculate the eigenvectors
    matrix::eigensystem_t es;
    m.SolveEigensystem(es);
    
    vec eigenv1=es.eigenvecs[0];
    vec eigenv2=es.eigenvecs[1];
    vec eigenv3=es.eigenvecs[2];
    
    _out->seteigenvec1(eigenv1);
    _out->seteigenvec2(eigenv2);
    _out->seteigenvec3(eigenv3);
    
    
    vec u = es.eigenvecs[0];
    vec v = _matrix[1]._in->getPos() - _matrix[0]._in->getPos();
    v.normalize();
    
    _out->setV(v);
    
    vec w = _matrix[2]._in->getPos() - _matrix[0]._in->getPos();
    w.normalize();
    
    //store the eigenvalues for the tensor of gyration
    double eigenvalue1 = es.eigenvalues[0];
    double eigenvalue2 = es.eigenvalues[1];
    double eigenvalue3 = es.eigenvalues[2];
    
    
    _out->setval1(eigenvalue1);
    _out->setval2(eigenvalue2);
    _out->setval3(eigenvalue3);
    
    if((v^w)*u < 0) u=vec(0.,0.,0.)-u;
    _out->setU(u);
    
    //write out w
    w=u^v;
    w.normalize();
    _out->setW(w);
    
    //out.BeadV(_out) = v;
    
    //out.BeadW(_out) = es.eigenvecs[2];
}
