// 
// File:   map.h
// Author: ruehle
//
// Created on April 13, 2007, 4:31 PM
//

#ifndef _map_H
#define	_map_H

#include <string>
#include <vector>
#include <tools/vec.h>
#include <tools/property.h>
#include "molecule.h"

using namespace std;


class BeadMap;
/*******************************************************
    Mapper class, collection of maps
*******************************************************/
class Map
{
public:
    Map(Molecule &in, Molecule &out)
        : _in(in), _out(out) {}
    ~Map();
    
    void AddBeadMap(BeadMap *bmap) { _maps.push_back(bmap); }

    void Apply();

protected:
    Molecule _in, _out;
    vector<BeadMap *> _maps;
};

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap
{
public:
    virtual ~BeadMap() {};
    virtual void Apply() = 0;
    virtual void Initialize(Molecule *in, Bead *out, Property *opts_map, Property *opts_bead);
protected:
    Molecule *_in;
    Bead *_out;
    Property *_opts_map;
    Property *_opts_bead;
};

inline void BeadMap::Initialize(Molecule *in, Bead *out, Property *opts_bead, Property *opts_map)
{
    _in = in; _out = out; _opts_map = opts_map; _opts_bead = opts_bead;
}

/*******************************************************
    Linear map for spherical beads
*******************************************************/
class Map_Sphere
    : public BeadMap
{
public:
    Map_Sphere() {}
    void Apply();

    void Initialize(Molecule *in, Bead *out, Property *opts_bead, Property *opts_map);

protected:
    void AddElem(Bead *in, double weight);

    struct element_t {
        Bead *_in;
        double _weight;
    };
    vector<element_t> _matrix;
};

inline void Map_Sphere::AddElem(Bead *in, double weight)
{
    element_t el;
    el._in = in;
    el._weight = weight;
    _matrix.push_back(el);
}

/*******************************************************
    Linear map for ellipsoidal bead
*******************************************************/
class Map_Ellipsoid
    : public Map_Sphere
{
public:
    Map_Ellipsoid() { }
    void Apply();
    
protected:
};


#endif	/* _map_H */

