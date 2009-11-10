// 
// File:   cgmoleculedef.h
// Author: ruehle
//
// Created on April 23, 2007, 5:58 PM
//

#ifndef _cgmoleculedef_H
#define	_cgmoleculedef_H

#include <string>
#include <vector>
#include <map>
#include <string>
#include <tools/property.h>
#include "map.h"
#include <tools/types.h>
#include "exclusionlist.h"
#include "molecule.h"

using namespace std;

/**
    \brief definition of a coarse grained molecule

    This class is to define a coarse grained molecule, which includes the topology, mapping, ...

    \todo clean up this class, do the bonded interactions right!!!!
    \todo check for consistency of xml file, seperate xml parser and class!!
*/
class CGMoleculeDef
{
public:
    CGMoleculeDef() {}
    ~CGMoleculeDef();
        
    Molecule *CreateMolecule(Topology & top);
    Map *CreateMap(Molecule &in, Molecule &out);

    void Load(string filename);
    
    const string &getName() { return _name; }
    const string &getIdent() { return _ident; }
    
    ExclusionList *CreateExclusionList(Molecule &atomistic);
    
private:
    Property _options;
    
    struct beaddef_t {
        string _name;
        string _type;
        byte_t _symmetry;
        string _mapping;
        vector<string> _subbeads;
        Property *_options;
    };    

    // name of the coarse grained molecule
    string _name;
    // name of the molecule to coarse grain
    string _ident;
    
    // beads of the cg molecule
    vector<beaddef_t *> _beads;
    map<string, beaddef_t *> _beads_by_name;
    
    // mapping schemes
    map<string, Property *> _maps;
    
    list<Property *> _bonded;
    
    void ParseTopology(Property &options);
    void ParseBeads(Property &options);
    void ParseBonded(Property &options);
    void ParseMapping(Property &options);
        
    beaddef_t *getBeadByName(const string &name);
    Property *getMapByName(const string &name);
};

#endif	/* _cgmoleculedef_H */

