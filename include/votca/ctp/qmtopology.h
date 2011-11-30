/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
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

namespace votca { namespace ctp {


/**
    \brief topology of conjugated segments and rigid fragments

    contains rigid fragments (centers of mass and orientations).
    Given a CG toplogy, it adds orientations to CG beads onto and
    updates positions of rigid fragments. Crg units should be associated to
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

    ///Initializes the charge units
    //void InitChargeUnits();
    
    /// \brief Create a new bead
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
    //void ComputeAllTransferIntegrals();

    /// find a conjugated segment by name
    QMCrgUnit *GetCrgUnitByName(const string &name);

    QMCrgUnit *getCrgUnit(int id);

    QMCrgUnit *CreateCrgUnit(const string &name, const string &type_name, int molid);
    QMCrgUnit *CreateCrgUnit(int id, const string &name, const string &type_name, int molid);

    /// find a charge unit type by name
    CrgUnitType *GetCrgUnitTypeByName(const string &type_name);

    //Copy charges to either charged or neutral case
    void CopyCharges(CrgUnit *crg, Molecule *mol);
    void CopyChargesOccupied(CrgUnit *crg, Molecule *mol);

    int getDatabaseId() { return _db_id; };
    void setDatabaseId(int id) { _db_id = id; }

protected:

    QMNBList _nblist;
    JCalc _jcalc;
    map <string, QMCrgUnit*> _mcharges;
    vector < QMCrgUnit *> _crgunits;
    map <int, QMCrgUnit *> _crgunits_by_id;

    int _db_id;
};

inline QMCrgUnit *QMTopology::getCrgUnit(int id)
{
    map<int, QMCrgUnit*>::iterator iter;
    iter = _crgunits_by_id.find(id);
    if(iter == _crgunits_by_id.end())
        throw std::runtime_error("did not find crgunit with id " + boost::lexical_cast<string>(id));
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

inline CrgUnitType *QMTopology::GetCrgUnitTypeByName(const string &type_name)
{
    return _jcalc.GetCrgUnitTypeByName(type_name);
}

}}

#endif	/* _CRGTOPOLOGY_H */

