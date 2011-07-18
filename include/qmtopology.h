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
    vector <QMCrgUnit *>& CrgUnits(){
        return _crgunits;
    }

    void Cleanup();
    /// update the topology based on cg positons
    void Update(Topology &cg_top);
    /// update the topology based on cg positons
    void Initialize(Topology &cg_top);
    ///Initialises the charge units
    void InitChargeUnits();
    
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

    QMCrgUnit *getCrgUnit(int id);

    QMCrgUnit *CreateCrgUnit(const string &name, const string &type_name, int molid);
    QMCrgUnit *CreateCrgUnit(int id, const string &name, const string &type_name, int molid);


    //Copy charges to either charged or neutral case
    void CopyCharges(CrgUnit *crg, Molecule *mol);
    void CopyChargesOccupied(CrgUnit *crg, Molecule *mol);

protected:

    QMNBList _nblist;
    JCalc _jcalc;
    map <string, QMCrgUnit*> _mcharges;
    vector < QMCrgUnit *> _crgunits;
    map <int, QMCrgUnit *> _crgunits_by_id;
};

inline QMCrgUnit *QMTopology::getCrgUnit(int id)
{
    map<int, QMCrgUnit*>::iterator iter;
    iter = _crgunits_by_id.find(id);
    if(iter == _crgunits_by_id.end())
        throw std::runtime_error("did not find crgunit with id " + lexical_cast<string>(id));
    return iter->second;
}


inline QMCrgUnit *QMTopology::GetCrgUnitByName(const string &name)
{
    map<string, QMCrgUnit *>::iterator iter;
    iter = _mcharges.find(name);
    if(iter!=_mcharges.end())
        return iter->second;
    return NULL;
}

#endif	/* _CRGTOPOLOGY_H */

