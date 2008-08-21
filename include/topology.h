// 
// File:   topology.h
// Author: ruehle
//
// Created on April 5, 2007, 11:35 AM
//

#ifndef _topology_H
#define	_topology_H

#include <vector>
#include <map>
#include <assert.h>
#include "beadinfo.h"
#include "moleculeinfo.h"
#include "residue.h"
#include "interaction.h"

typedef vector<MoleculeInfo *> MoleculeContainer;
typedef vector<BeadInfo *> BeadContainer;
typedef vector<Residue *> ResidueContainer;
typedef vector<Interaction *> InteractionContainer;


using namespace std;

/**
    \brief topology of a molecule

    The CTopology class stores the topology of the system like the beads, bonds, molecules and residues.

    \todo internal management for ids and indices
 **/
class Topology
{
public:
    /// constructor
    Topology() {}
    
    ~Topology();
    
    /// cleans up all the stored data
    void Cleanup();
    
    /// creates a new Bead
    BeadInfo *CreateBead(byte_t symmetry, string name, int resnr, double m, double q);
    //int AddBead(CBeadInfo *bead);

    /// creates a new molecule
    MoleculeInfo *CreateMolecule(string name);
    
    /// creates a new residue
    Residue *CreateResidue(string name);
    
    /** 
        \brief create molecules based on the residue
        This function scans the topology and creates molecules based on the resiude id.
        All beads with the same resid are put int one molecule.
    */
    void CreateMoleculesByResidue();
    
    /** 
        \brief put the whole topology in one molecule
        This function creates one big molecule for all beads in the topology.
    */
    void CreateOneBigMolecule(string name);
    
    /// number of molecules in the system
    int MoleculeCount() { return _molecules.size(); }
    /// number of beads in the system
    int BeadCount() { return _beads.size(); }
    
    /// number of beads in the system
    int ResidueCount() { return _residues.size(); }
       
    /// get a molecule by an index
    MoleculeInfo *MoleculeByIndex(int index);
    
    /// get the bead container
    BeadContainer &getBeads() { return _beads; }
    /// get the molecule container
    MoleculeContainer &getMolecules() { return _molecules; }

    /// get the bonded interaction container
    InteractionContainer &getBondedInteractions() { return _interactions; }
    
    
    // \todo change AddBondedinteraction to Create bonded interaction, that only topology can create interactions
    void AddBondedInteraction(Interaction *ic);
    
    BeadInfo *getBead(const int i) { return _beads[i]; }
    Residue *getResidue(const int i) { return _residues[i]; }
    MoleculeInfo *getMolecule(const int i) { return _molecules[i]; }
   
    
    /// a hack to get vale's topology to read
    void ClearMoleculeList(){
        _molecules.clear();
    }
    
    void Add(Topology *top);
    
    void RenameMolecules(string range, string name);
    
private:
    /// beads in the topology
    BeadContainer _beads;
    
    /// molecules in the topology
    MoleculeContainer _molecules;
    
    /// residues in the topology
    ResidueContainer _residues;
    
    /// bonded interactions in the topology
    InteractionContainer _interactions;
    
    map<string, int> _interaction_groups;
};

inline BeadInfo *Topology::CreateBead(byte_t symmetry, string name, int resnr, double m, double q)
{
    
    BeadInfo *b = new BeadInfo(_beads.size(), symmetry, name, resnr, m, q);    
    _beads.push_back(b);
    return b;
}

inline MoleculeInfo *Topology::CreateMolecule(string name)
{
    MoleculeInfo *mol = new MoleculeInfo(_molecules.size(), name);
    _molecules.push_back(mol);
    return mol;
}

inline Residue *Topology::CreateResidue(string name)
{
    Residue *res = new Residue(_molecules.size(), name);
    _residues.push_back(res);
    return res;
}

inline MoleculeInfo *Topology::MoleculeByIndex(int index)
{
    return _molecules[index];
}

#endif	/* _topology_H */

