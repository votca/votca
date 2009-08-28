// 
// File:   cgmoleculedef.cc
// Author: ruehle
//
// Created on April 23, 2007, 6:36 PM
//

#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include "topology.h"
#include "cgmoleculedef.h"
#include <tools/tokenizer.h> 
#include "interaction.h"

using boost::lexical_cast;

CGMoleculeDef::~CGMoleculeDef()
{
    {
        vector<beaddef_t *>::iterator i;
        for(i=_beads.begin();i!=_beads.end();++i)
            delete *i;
        _beads.clear();
    }
    {
        map<string, mapdef_t *>::iterator i;
        for(i=_maps.begin();i!=_maps.end();++i)
            delete (*i).second;
        _maps.clear();
    }
    {
        vector<forcedef_t *>::iterator i;
        for(i=_bonded.begin();i!=_bonded.end();++i)
            delete *i;
        _bonded.clear();
    }
}

void CGMoleculeDef::Load(string filename)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;
    
    // open the xml file
    //cout << filename << endl;
    doc = xmlParseFile(filename.c_str());
    if(doc == NULL) 
        throw runtime_error("Error on open xml bead map: " + filename);
    
    node = xmlDocGetRootElement(doc);
    
    if(node == NULL) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, empty xml document: " + filename);
    }
    
    if(xmlStrcmp(node->name, (const xmlChar *)"cg_molecule")) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, in xml file: " + filename);
    }
    
    // parse xml tree
    for(node = node->xmlChildrenNode; node != NULL; node = node->next) {
            key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
            if(!xmlStrcmp(node->name, (const xmlChar *) "name"))
                _name = reinterpret_cast<char *>(key);
            else if(!xmlStrcmp(node->name, (const xmlChar *) "ident"))
                _ident = reinterpret_cast<char *>(key);
            else if(!xmlStrcmp(node->name, (const xmlChar *) "topology"))
                ParseTopology(doc, node->xmlChildrenNode);
            else if(!xmlStrcmp(node->name, (const xmlChar *) "maps"))
                ParseMapping(doc, node->xmlChildrenNode);
            xmlFree(key);
    }   
    
    // free the document
    //xmlFree(node);
    xmlFreeDoc(doc);
    xmlCleanupParser();
}

void CGMoleculeDef::ParseTopology(xmlDocPtr doc, xmlNodePtr cur)
{
    for(; cur != NULL; cur = cur->next) {
        if(!xmlStrcmp(cur->name, (const xmlChar *)"cg_beads"))
            ParseBeads(doc, cur->xmlChildrenNode);
        if(!xmlStrcmp(cur->name, (const xmlChar *)"cg_bonded"))
            ParseBonded(doc, cur->xmlChildrenNode);
    }
}

void CGMoleculeDef::ParseBeads(xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;
    xmlNodePtr node;
    
    for(; cur != NULL; cur = cur->next) {
        if(!xmlStrcmp(cur->name, (const xmlChar *)"cg_bead")) {
            beaddef_t *beaddef = new beaddef_t;
            beaddef->_symmetry = 1;
            for(node = cur->xmlChildrenNode; node != NULL; node = node->next) {
                key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
                if(!xmlStrcmp(node->name, (const xmlChar *) "name"))
                    beaddef->_name = reinterpret_cast<char *>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "type"))
                    beaddef->_type = reinterpret_cast<char *>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "symmetry"))
                    beaddef->_symmetry = lexical_cast<int>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "mapping"))
                    beaddef->_mapping = reinterpret_cast<char *>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "beads")) {
                    string s(reinterpret_cast<char *>(key));
                    Tokenizer tok(s, " \n\t");                   
                    for(Tokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg)
                        beaddef->_subbeads.push_back(*beg);
                }
                else if (node->type == XML_ELEMENT_NODE) {
                    ParseNode(beaddef->_misc[(const char *)node->name], doc, node);
                }
                xmlFree(key);
            }            
            _beads.push_back(beaddef);
            _beads_by_name[beaddef->_name] = beaddef;
        }
    }
}

void CGMoleculeDef::ParseBonded(xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;
    xmlNodePtr node;
    
    for(; cur != NULL; cur = cur->next) {
        if( (!xmlStrcmp(cur->name, (const xmlChar *)"bond")) 
            || (!xmlStrcmp(cur->name, (const xmlChar *)"angle")) || (!xmlStrcmp(cur->name, (const xmlChar *)"dihedral")) ) {
            forcedef_t *fdef = new forcedef_t;            
            fdef->_type = reinterpret_cast<const char *>(cur->name); //"bond";
            for(node = cur->xmlChildrenNode; node != NULL; node = node->next) {
                key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
                if(!xmlStrcmp(node->name, (const xmlChar *) "name"))
                    fdef->_name = reinterpret_cast<char *>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "beads")) {
                    string s(reinterpret_cast<char *>(key));
                    Tokenizer tok(s, " \n\t");
                    for(Tokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg) {
                        //if((!xmlStrcmp(cur->name, (const xmlChar *)"dihedral"))) cout << *beg << endl;
                        fdef->_atoms.push_back(*beg);
                    }
                }
                else if (node->type == XML_ELEMENT_NODE) {
                    ParseNode(fdef->_misc[(const char *)node->name], doc, node);
                }
                xmlFree(key);
            }            
            _bonded.push_back(fdef);
        }        
    }
}

void CGMoleculeDef::ParseMapping(xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;
    xmlNodePtr node;

    for(; cur != NULL; cur = cur->next) {
        if(!xmlStrcmp(cur->name, (const xmlChar *)"map")) {
            mapdef_t *mapdef = new mapdef_t;
            
            for(node = cur->xmlChildrenNode; node != NULL; node = node->next) {
                key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
                if(!xmlStrcmp(node->name, (const xmlChar *) "name"))
                    mapdef->_name = reinterpret_cast<char *>(key);
                else if(!xmlStrcmp(node->name, (const xmlChar *) "weights")) {
                    double m=0;
                    string s(reinterpret_cast<char *>(key));
                    Tokenizer tok(s, " \n\t");
                    for(Tokenizer::iterator beg=tok.begin(); beg!=tok.end();++beg) {
                        double tmp = lexical_cast<double>(*beg);
                        m+=tmp;
                        mapdef->_weights.push_back(tmp);
                    }
                    for(size_t i=0; i<mapdef->_weights.size(); i++)
                        mapdef->_weights[i]/=m;
                }
                else if (node->type == XML_ELEMENT_NODE) {
                    ParseNode(mapdef->_misc[(const char *)node->name], doc, node);
                }
                xmlFree(key);
            }
            _maps[mapdef->_name] = mapdef;            
        }
    }
}

void CGMoleculeDef::ParseNode(option_t &op, xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;
    xmlNodePtr node;

    op._name = (const char *)cur->name;
    key = xmlNodeListGetString(doc, cur, 1);
    op._value = (const char *)key;
    xmlFree(key);
    
    for(node = cur->xmlChildrenNode; node != NULL; node = node->next) {
        if (node->type == XML_ELEMENT_NODE)
            ParseNode(op._childs[(const char *)node->name], doc, node->xmlChildrenNode);            
    }
}

Molecule * CGMoleculeDef::CreateMolecule(Topology & top)
{   
    // add the residue names
    Residue *res = top.CreateResidue(_name);
    Molecule *minfo = top.CreateMolecule(_name);
    
    // create the atoms
    vector<beaddef_t *>::iterator iter;
    for(iter = _beads.begin(); iter != _beads.end(); ++iter) {
        Bead *bead;
        BeadType *bt = top.GetOrCreateBeadType((*iter)->_type);
        bead = top.CreateBead((*iter)->_symmetry, (*iter)->_name, bt, res->getId(), 0, 0);
        minfo->AddBead(bead, bead->getName());        
        
    }    
    
    // create the bonds
    vector<forcedef_t *>::iterator fdef;
    for(fdef = _bonded.begin(); fdef!=_bonded.end(); ++fdef) {
        if((*fdef)->_type == "bond") {
            for(size_t i=0; i<(*fdef)->_atoms.size(); i+=2) {
                int atom1, atom2;
                atom1 = minfo->getBeadIdByName((*fdef)->_atoms[i]);
                atom2 = minfo->getBeadIdByName((*fdef)->_atoms[i+1]);
                if(atom1 < 0)
                    runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i] + " not found"));
                if(atom2 < 0) 
                    throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+1] + " not found"));
                IBond *ib = new IBond(atom1, atom2);
                ib->setGroup((*fdef)->_name);
                ib->setIndex(i/2);
                ib->setMolecule(minfo->getId());
                top.AddBondedInteraction(dynamic_cast<Interaction *>(ib));
            }
        }
        else if((*fdef)->_type == "angle") {
            for(size_t i=0; i<(*fdef)->_atoms.size(); i+=3) {
                int atom1, atom2, atom3;
                atom1 = minfo->getBeadIdByName((*fdef)->_atoms[i]);
                atom2 = minfo->getBeadIdByName((*fdef)->_atoms[i+1]);
                atom3 = minfo->getBeadIdByName((*fdef)->_atoms[i+2]);
                
                if(atom1 < 0)
                    throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i] + " not found"));
                if(atom2 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+1] + " not found"));
                if(atom3 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+3] + " not found"));
                IAngle *ib = new IAngle(atom1, atom2, atom3);
                ib->setGroup((*fdef)->_name);
                ib->setIndex(i/3);
                ib->setMolecule(minfo->getId());
                top.AddBondedInteraction(dynamic_cast<Interaction *>(ib));
            }
        }
        else if((*fdef)->_type == "dihedral") {
            for(size_t i=0; i<(*fdef)->_atoms.size(); i+=4) {
                int atom1, atom2, atom3, atom4;
                atom1 = minfo->getBeadIdByName((*fdef)->_atoms[i]);
                atom2 = minfo->getBeadIdByName((*fdef)->_atoms[i+1]);
                atom3 = minfo->getBeadIdByName((*fdef)->_atoms[i+2]);
                atom4 = minfo->getBeadIdByName((*fdef)->_atoms[i+3]);
                if(atom1 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i] + " not found"));
                if(atom2 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+1] + " not found"));
                if(atom3 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+3] + " not found"));
                if(atom4 < 0) throw runtime_error(string("error while trying to create bond, bead " + (*fdef)->_atoms[i+3] + " not found"));
                IDihedral *ib = new IDihedral(atom1, atom2, atom3, atom4);
                ib->setGroup((*fdef)->_name);
                ib->setIndex(i/4);
                ib->setMolecule(minfo->getId());
                top.AddBondedInteraction(dynamic_cast<Interaction *>(ib));
            }
        }
    }
    return minfo;
}

/// \bug crashes when map does not exist!
Map *CGMoleculeDef::CreateMap(Molecule &in, Molecule &out)
{
    vector<beaddef_t *>::iterator def;
    mapdef_t *mdef;
    Map *map;
    int iout, iin;
    /// \bug check weather molecule and #define match
    map = new Map();
    for(def = _beads.begin(); def != _beads.end(); ++def) {
        vector<int>::iterator iter;
        iout = out.getBeadByName((*def)->_name);
        if(iout < 0) 
            throw runtime_error(string("mapping error: molecule " + (*def)->_name + " does not exist"));
        
        mdef = getMapByName((*def)->_mapping);
        if(!mdef)
            throw runtime_error(string("mapping " + (*def)->_mapping + " not found"));
        
        if((*def)->_subbeads.size() != mdef->_weights.size())
            throw runtime_error(string("number of subbeads in " + (*def)->_name + "and number of weights in map " + (*def)->_mapping +" do not match"));
        
        switch((*def)->_symmetry) {
        case 1:
        {
            Map_Sphere *bmap = new Map_Sphere(iout);           
            for(size_t i=0; i < (*def)->_subbeads.size(); ++i) {
                iin = in.getBeadByName((*def)->_subbeads[i]);
                if(iin < 0) 
                    throw runtime_error(string("mapping error: molecule " + (*def)->_subbeads[i] + " does not exist"));
                bmap->AddElem(iin, mdef->_weights[i]);
            }
            map->AddBeadMap(bmap);
        }
            break;        
        case 3:
        {
            Map_Ellipsoid *bmap = new Map_Ellipsoid(iout); ;
            for(size_t i=0; i < (*def)->_subbeads.size(); ++i) {
                iin = in.getBeadByName((*def)->_subbeads[i]);
                if(iin < 0) 
                    throw runtime_error(string("mapping error: molecule " + (*def)->_subbeads[i] + " does not exist"));
                bmap->AddElem(iin, mdef->_weights[i]);
            }
            map->AddBeadMap(bmap);
        }
            break;
        default:
            throw runtime_error(string("unknown symmetry in bead definition!"));
        }
    }
    return map;
}

/**
 * \todo Check this function for multiple molecules!!!!!!! 
 */
ExclusionList *CGMoleculeDef::CreateExclusionList(Molecule &atomistic)
{
    list<int> exclude;
    int natoms;
    
    ExclusionList *ex = new ExclusionList();
    ex->ExcludeAll(atomistic.BeadCount());        
    vector<forcedef_t *>::iterator iter;
    
    // reintroduce bead internal nonbonded interaction
    vector<beaddef_t *>::iterator bd_iter;
    for(bd_iter = _beads.begin(); bd_iter!=_beads.end(); ++bd_iter) {
        exclude.clear();
        beaddef_t *bead = *bd_iter;
        for(vector<string>::iterator sb=bead->_subbeads.begin(); sb!=bead->_subbeads.end(); ++sb) {
                    exclude.push_back(atomistic.getBeadId(atomistic.getBeadIdByName(*sb)));
//                    cout << atomistic.getBeadId(atomistic.getBeadByName(*sb)) << " ";
        }
//        cout << endl;
        ex->Remove(exclude);
    }
    
    // reintroduce nonbonded interactions for bonded beads
    for(iter = _bonded.begin(); iter!=_bonded.end(); ++iter) {
        if((*iter)->_type == "bond")
            natoms = 2;              
        else if((*iter)->_type == "angle")
            natoms = 3;              
        else if((*iter)->_type == "dihedral") {
            natoms = 4;              
        }
        else throw runtime_error(string("unknown bond type"));
        
        for(size_t i=0; i<(*iter)->_atoms.size(); i+=natoms) {
            exclude.clear();
        
            for(int j=0; j<natoms; j++) {
                beaddef_t *bead = getBeadByName((*iter)->_atoms[i+j]);
                if(bead == NULL) 
                    throw runtime_error(string("error while trying to create exclusion list, bead " + (*iter)->_atoms[i+j] + " not found"));
                for(vector<string>::iterator sb=bead->_subbeads.begin(); sb!=bead->_subbeads.end(); ++sb) {
                    exclude.push_back(atomistic.getBeadId(atomistic.getBeadIdByName(*sb)));
                }                
            }
            ex->Remove(exclude);
        }
    }
   
    return ex;
}

CGMoleculeDef::beaddef_t *CGMoleculeDef::getBeadByName(const string &name)
{
    map<string, beaddef_t*>::iterator iter = _beads_by_name.find(name);
    if(iter == _beads_by_name.end()) {
        std::cout << "cannot find: <" << name << "> in " << _name << "\n";
        return NULL;        
    }
    //assert(iter != _beadmap.end());
    //return (*iter).second;
    return (*iter).second;
}

CGMoleculeDef::mapdef_t *CGMoleculeDef::getMapByName(const string &name)
{
    map<string, mapdef_t*>::iterator iter = _maps.find(name);
    if(iter == _maps.end()) {
        std::cout << "cannot find map " << name << "\n";
        return NULL;        
    }
    //assert(iter != _beadmap.end());
    //return (*iter).second;
    return (*iter).second;
}

void CGMoleculeDef::OutputNode(option_t &opt)
{
    map<string,option_t>::iterator i;
    cout << "<<\n" << opt._name << ": " << opt._value << "\n";
    for(i=opt._childs.begin(); i!=opt._childs.end();++i)
        OutputNode((*i).second);
    cout << ">>\n";
}
    
