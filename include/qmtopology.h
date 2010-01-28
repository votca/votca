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
#include <moo/crgunittype.h>
#include <moo/crgunit.h>
#include <moo/jcalc.h>
#include <moo/units.h>

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
    list < CrgUnit *>& crglist(){
        return _lcharges;
    }
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
    void AddAtomisticBeads(CrgUnit * crg, Topology * totop);

    /// computes all transfer integrals (note that the nblist must be initialised by the user!)
    void ComputeAllTransferIntegrals();

    /// comutes all the electrostatic energies
    void ComputeAllElectrostaticEnergies(const double &epsilon=3.5);

    /// find a crg unit by name
    CrgUnit *GetCrgUnitByName(const string &name);

    CrgUnit *CreateCrgUnit(const string &name, const string &type_name, int molid);
    
protected:

    QMNBList _nblist;
    JCalc _jcalc;
    map <string, CrgUnit*> _mcharges;
    list < CrgUnit *> _lcharges;
    
    ///Initialises the charge units
    void InitChargeUnits();
};

inline CrgUnit *QMTopology::GetCrgUnitByName(const string &name)
{
    map<string, CrgUnit *>::iterator iter;
    iter = _mcharges.find(name);
    if(iter!=_mcharges.end())
        return iter->second;
    return NULL;
}

inline CrgUnit *QMTopology::CreateCrgUnit(const string &name, const string &type_name, int molid)
{
    if(GetCrgUnitByName(name))
        throw std::runtime_error("charge unit with name " + name + " already exists");
    CrgUnit *crg;
    crg = _jcalc.CreateCrgUnit(_lcharges.size(), type_name, molid);
    _mcharges.insert(make_pair(name, crg));
    _lcharges.push_back(crg);
    crg->setName(name);
    return crg;
}

#endif	/* _CRGTOPOLOGY_H */

