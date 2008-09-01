/// \addtogroup csg
///@{
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
#include "molecule.h"

using namespace std;


class BeadMap;
/*******************************************************
    Mapper class, collection of maps
*******************************************************/
class Map
{
public:
    Map() {}
    ~Map();
    void Apply(Molecule &in, Molecule &out);
    
    void AddBeadMap(BeadMap *bmap) { _maps.push_back(bmap); }
protected:
    vector<BeadMap *> _maps;
};

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap
{
public:
    virtual ~BeadMap() {};
    virtual void Apply(Molecule &in, Molecule &out) = 0;
    
    /*
        \brief clone map object
        Clones a map object. Each mapper class must implement this function to create an instance of itself.
        It is needed by the mapper factory to create a map of the given type.
    */
    //virtual BeadMap *Clone() = 0;
        
};

/*******************************************************
    Linear map for spherical beads
*******************************************************/
class Map_Sphere
    : public BeadMap
{
public:
    Map_Sphere(int out) : _out(out) {};
    void Apply(Molecule &in, Molecule &out);

    void AddElem(int in, double weight);
    
    //BeadMap *Clone() { return new Map_Sphere(_out); }
    
protected:
    Map_Sphere() {}
    
    struct element_t {
        int _in;
        double _weight;
    };
    int _out;
    vector<element_t> _matrix;
};

inline void Map_Sphere::AddElem(int in, double weight)
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
    Map_Ellipsoid(int out) : _out(out) { }
    void Apply(Molecule &in, Molecule &out);

    //void AddElem(int out, int in, double weight);
    //BeadMap *Clone() { return new Map_Ellipsoid(_out); }
    
protected:
    int _out;
    /*struct element_t {
        int _in, _out;
        double _weight;
    };
    
    vector<element_t> _matrix;*/
};


#endif	/* _map_H */

/// @}
