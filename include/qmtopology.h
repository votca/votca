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

class QMNBList;

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

    //assign neighbours
    void AssignNBList(QMNBList * nblist);

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
protected:

    QMNBList *_nblist;
    JCalc _jcalc;
    map <string, CrgUnit*> _mcharges;
    list < CrgUnit *> _lcharges;
};

inline  void QMTopology::AssignNBList(QMNBList  * nblist){
    _nblist= nblist;
}


#endif	/* _CRGTOPOLOGY_H */

