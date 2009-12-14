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

#include "qmbead.h"

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

    ///at each evaluate CG step we will need to reassess the QMBeads
    //int UpdateQMTopology();

    /// update the topology based on cg positons
    void Update(Topology &cg_top);

    /// \brief Cretae a new bead
    /// We overload CreateBead to create QMBead, this is needed to make
    /// CopyTopologyData work
    Bead *CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q);

    void LoadListCharges(const string &file);
protected:

    NBList *_nblist;
    JCalc _jcalc;

};

inline Bead *QMTopology::CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q)
{
    QMBead *b = new QMBead(this, _beads.size(), type, symmetry, name, resnr, m, q);
    _beads.push_back(b);
    //initialise the crgunit *
    // NOTE: I cannot find this famours getOption in bead?!
    string namecrgunittype = bead->getType()->getName();
    string intpos = bead->getOptions->get("qm.position".as<int> ());
    string namecrgunit = bead->getOptions->get("qm.crgunittype".as<string> ());
    
    CrgUnitType * _crtgtype  = _jcalc.GetCrgUnitTypeByName(namecrgunittype);

    //determine whether it  has been created already
    /does the combination of bead->Molecule()->molID + namecrgunit exist?
        yes-> do nothing
        no -> create a crgunit of tupe crgtype, with molid bla bla and pos xyz
    return b;
}

#endif	/* _CRGTOPOLOGY_H */

