/* 
 * File:   crgtopology.h
 * Author: james
 *
 * Created on 12 December 2009, 13:50
 */

#ifndef _QMTOPOLOGY_H
#define	_QMTOPOLOGY_H

#include <votca/csg/topology.h>
#include <votca/csg/nblist.h>
#include <votca/moo/crgunittype.h>
#include "qmcrgunit.h"
#include <votca/moo/jcalc.h>
#include <votca/moo/units.h>

#include "qmbead.h"
#include "qmnblist.h"


/**
    \brief topology of qmbeads

    contains the beads describing the c.o.m. of each cahrge unit
    given a CG toplogy, it should be apply the mapping for the cg beads -> qm beads and
    it should update the position of the crg unit. Crg units should be associated to
    these qm beads and not to any other.
*/


class QMTopology : public Topology
{
public:
    QMTopology();
    ~QMTopology();

    /// QM neighbor list [nm]
    QMNBList & nblist(){
        return _nblist;
    }
    /// List of charge units [Bohr]
    vector < QMCrgUnit *>& CrgUnits(){
        return _crgunits;
    }

    void Cleanup();
    /// update the topology based on cg positons
    void Initialize(Topology &cg_top);
    /// update the topology based on cg positons
    void Update(Topology &cg_top);

    /// \brief Cretae a new bead
    /// We overload CreateBead to create QMBead, this is needed to make
    /// CopyTopologyData work
    Bead *CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q);

    /// Loads the definitions of crgunits into jcalc
    void LoadListCharges(const string &file);

    /// Returns the underlying crgunit type lists etc
    JCalc &GetJCalc() { return _jcalc; }

    ///Loads the atomistic beads (from mol_and_orb) into totop from the CrgUnit defined by namecrgunit and molid
    Molecule *AddAtomisticBeads(CrgUnit * crg, Topology * totop);

    /// computes all transfer integrals (note that the nblist must be initialised by the user!)
    // TODO: this function should not be in qmtopology!
    void ComputeAllTransferIntegrals();

    /// find a crg unit by name
    QMCrgUnit *GetCrgUnitByName(const string &name);

    QMCrgUnit *GetCrgUnit(int index);

    QMCrgUnit *CreateCrgUnit(const string &name, const string &type_name, int molid);


    //Copy charges to either charged or neutral case
    void CopyCharges(CrgUnit *crg, Molecule *mol);
    void CopyChargesOccupied(CrgUnit *crg, Molecule *mol);

protected:

    QMNBList _nblist;
    JCalc _jcalc;
    map <string, QMCrgUnit*> _mcharges;
    vector < QMCrgUnit *> _crgunits;
    
    ///Initialises the charge units
    void InitChargeUnits();
};

inline QMCrgUnit *QMTopology::GetCrgUnit(int index)
{
    if(index >= _crgunits.size())
        throw std::runtime_error("error, crgunit index out of bounds");
    return _crgunits[index];
}


inline QMCrgUnit *QMTopology::GetCrgUnitByName(const string &name)
{
    map<string, QMCrgUnit *>::iterator iter;
    iter = _mcharges.find(name);
    if(iter!=_mcharges.end())
        return iter->second;
    return NULL;
}

inline QMCrgUnit *QMTopology::CreateCrgUnit(const string &name, const string &type_name, int molid)
{
    if(GetCrgUnitByName(name))
        throw std::runtime_error("charge unit with name " + name + " already exists");
    QMCrgUnit *crg;

    CrgUnitType *type = _jcalc.GetCrgUnitTypeByName(type_name);
    if(!type)
        throw runtime_error("Charge unit type not found: " + type_name);
       
    crg = new QMCrgUnit(_crgunits.size(), type, molid);

    _mcharges.insert(make_pair(name, crg));
    _crgunits.push_back(crg);
    crg->setName(name);
    return crg;
}

#endif	/* _CRGTOPOLOGY_H */

