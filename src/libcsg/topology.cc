// 
// File:   topology.cc
// Author: ruehle
//
// Created on April 5, 2007, 12:30 PM
//

#include "topology.h"
#include <tools/rangeparser.h>

Topology::~Topology()
{
    Cleanup();
}

void Topology::Cleanup()
{
    // cleanup beads
    {
        BeadContainer::iterator i;
        for(i=_beads.begin();i<_beads.end();++i)
            delete *i;
        _beads.clear();
    }
    // cleanup molecules
    {
        MoleculeContainer::iterator i;
        for(i=_molecules.begin();i<_molecules.end();++i)
            delete *i;
        _molecules.clear();
    }
    // cleanup residues
    {
        ResidueContainer::iterator i;
        for(i=_residues.begin();i<_residues.end();++i)
            delete (*i);
        _residues.clear();
    }
    // cleanup interactions
    {
        InteractionContainer::iterator i;
        for(i=_interactions.begin();i<_interactions.end();++i)
            delete (*i);
        _interactions.clear();
    }
}

/// \todo clean up CreateMoleculesByResidue!
void Topology::CreateMoleculesByResidue()
{
    // first create a molecule for each residue    
    ResidueContainer::iterator res;    
    for(res=_residues.begin(); res!=_residues.end(); ++res) {
        CreateMolecule((*res)->getName());
    }
    
    // add the beads to the corresponding molecules based on their resid
    BeadContainer::iterator bead;    
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
        //MoleculeByIndex((*bead)->getResnr())->AddBead((*bead)->getId(), (*bead)->getName());

        MoleculeByIndex((*bead)->getResnr())->AddBead((*bead)->getId(), string("1:TRI:") + (*bead)->getName());
    }
    
    /// \todo sort beads in molecules that all beads are stored in the same order. This is needed for the mapping!
}

void Topology::CreateOneBigMolecule(string name)
{
    MoleculeInfo *mi = CreateMolecule(name);
    
    BeadContainer::iterator bead;    
    
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
        stringstream n("");
        n << (*bead)->getResnr() +1 << ":" <<  _residues[(*bead)->getResnr()]->getName() << ":" << (*bead)->getName();
        //cout << n.str() << endl;
        mi->AddBead((*bead)->getId(), n.str());
    }    
}

void Topology::Add(Topology *top)
{
    BeadContainer::iterator bead;
    ResidueContainer::iterator res;
    MoleculeContainer::iterator mol;
    
    int bead0=BeadCount();
    int res0=ResidueCount();
    int mol0=MoleculeCount();
    
    for(bead=top->_beads.begin(); bead!=top->_beads.end(); ++bead) {
        BeadInfo *bi = *bead;
        CreateBead(bi->getSymmetry(), bi->getName(), 0, bi->getResnr()+res0, bi->getM(), bi->getQ());
    }
    
    for(res=top->_residues.begin();res!=top->_residues.end(); ++res) {
        CreateResidue((*res)->getName());
    }
  
    // \todo beadnames in molecules!!
    for(mol=top->_molecules.begin();mol!=top->_molecules.end(); ++mol) {
        MoleculeInfo *mi = CreateMolecule((*mol)->getName());
        for(int i=0; i<mi->BeadCount(); i++) {
            mi->AddBead(mi->getBeadId(i), "invalid");
        }
    }
}

void Topology::RenameMolecules(string range, string name)
{
    RangeParser rp;
    RangeParser::iterator i;
    
    rp.Parse(range);
    for(i=rp.begin();i!=rp.end();++i) {
        getMolecule(*i-1)->setName(name);
    }
}

void Topology::AddBondedInteraction(Interaction *ic)
{
    map<string,int>::iterator iter;
    iter = _interaction_groups.find(ic->getGroup());
    if(iter!=_interaction_groups.end())
        ic->setGroupId((*iter).second);
    else {
        int i= _interaction_groups.size();
        _interaction_groups[ic->getGroup()] = i;
        ic->setGroupId(i);
    }
    _interactions.push_back(ic);
}


