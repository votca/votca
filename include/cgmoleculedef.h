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
#include <libxml/parser.h>
#include "map.h"
#include <tools/types.h>
#include "exclusionlist.h"
#include "moleculeinfo.h"

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
        
    MoleculeInfo *CreateMolecule(Topology & top);
    Map *CreateMap(MoleculeInfo &in, MoleculeInfo &out);

    void Load(string filename);
    
    const string &getName() { return _name; }
    const string &getIdent() { return _ident; }
    
    ExclusionList *CreateExclusionList(MoleculeInfo &atomistic);
    
    struct option_t {
       string _name;
       string _value;
       map<string,option_t> _childs;       
    };
    
private:
    struct beaddef_t {
        string _name;
        string _type;
        byte_t _symmetry;
        string _mapping;
        vector<string> _subbeads;
        map<string,option_t> _misc;
    };
    
    struct mapdef_t {
        string _name;
        vector<double> _weights;
        map<string,option_t> _misc;
    };

    struct forcedef_t {
        string _type;
        string _name;
        vector< string > _atoms;
        map<string,option_t> _misc;
    };

    // name of the coarse grained molecule
    string _name;
    // name of the molecule to coarse grain
    string _ident;
    
    // beads of the cg molecule
    vector<beaddef_t *> _beads;
    map<string, beaddef_t *> _beads_by_name;
    
    // mapping schemes
    map<string, mapdef_t *> _maps;
    
    vector<forcedef_t *> _bonded;
    
    void ParseTopology(xmlDocPtr doc, xmlNodePtr cur);
    void ParseBeads(xmlDocPtr doc, xmlNodePtr cur);
    void ParseBonded(xmlDocPtr doc, xmlNodePtr cur);
    void ParseMapping(xmlDocPtr doc, xmlNodePtr cur);
    
    void ParseNode(option_t &opt, xmlDocPtr doc, xmlNodePtr cur);
    void OutputNode(option_t &opt);
    
    beaddef_t *getBeadByName(const string &name);
    mapdef_t *getMapByName(const string &name);    
};

#endif	/* _cgmoleculedef_H */

