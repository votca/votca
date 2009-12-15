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
    void Initialize(Topology &cg_top);
    /// update the topology based on cg positons
    void Update(Topology &cg_top);

    /// \brief Cretae a new bead
    /// We overload CreateBead to create QMBead, this is needed to make
    /// CopyTopologyData work
    Bead *CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q);

    void LoadListCharges(const string &file);

    JCalc &GetJCalc() { return _jcalc; }
protected:

    NBList *_nblist;
    JCalc _jcalc;
    map <string, CrgUnit*> _mcharges;
    list < CrgUnit *> _lcharges;
};

inline Bead *QMTopology::CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q)
{
    QMBead *bead = new QMBead(this, _beads.size(), type, symmetry, name, resnr, m, q);
    _beads.push_back(bead);

    //initialise the crgunit * only if appropriate extra info is in the cg.xml file
    if (! (bead->Options()).exists("qm.crgunitname")){
        string namecrgunittype = bead->getType()->getName();
        int intpos = (bead->Options()).get("qm.position").as<int>();
        string namecrgunit = (bead->Options()).get("qm.crgunitname").as<string>();

        CrgUnitType *crgtype  = _jcalc.GetCrgUnitTypeByName(namecrgunittype);

        //determine whether it  has been created already
        int molid= bead->getMolecule()->getId();
        string molandtype = lexical_cast<string>(molid)+":"+namecrgunit;
        map <string, CrgUnit*>::iterator  itm= _mcharges.find(molandtype);
        if (itm != _mcharges.end()){
            vector <vec> empty;
            CrgUnit * acrg = new CrgUnit(empty, empty, empty, // this is because i dont want to cannot init all values at once
                _lcharges.size(), crgtype, molid);
            _mcharges.insert(make_pair(molandtype, acrg));
            _lcharges.push_back(acrg);
            bead->SetCrg(acrg);
            bead->SetiPos(intpos);
        }
        else{
            bead->SetCrg(NULL);
        }
    }
    return bead;
}

#endif	/* _CRGTOPOLOGY_H */

