#ifndef _QMNBList_H
#define	_QMNBList_H

#include "qmpair.h"
#include <votca/csg/pairlist.h>
#include <votca/csg/beadlist.h>
#include <moo/crgunit.h>

using namespace votca::csg;

class QMTopology;
using namespace votca::tools;
using namespace votca::csg;

/**
 * \brief Neighbour search for crg units
 *
 * This class wraps the NBList from csg to work on CrgUnits,
 * this all looks a bit cumbersome now, but will make things
 * nice if one want's to switch in between nbsearch algorithms
 *
 * */

class QMNBList
    : public PairList<QMCrgUnit *, QMPair>
{
public:
     QMNBList() {};
     ~QMNBList(){
         PairList<QMCrgUnit *, QMPair>::Cleanup();
     }
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }

    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() {return _cutoff; }

protected:

    bool Match(Bead *b1, Bead *b2, const vec &r, const double notused);
    double _cutoff;
    QMTopology *_father;

};

#endif	/* _QMNBList_H */

