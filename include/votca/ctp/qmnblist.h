#ifndef _QMNBList_H
#define	_QMNBList_H


#include <stdlib.h>
#include <votca/csg/pairlist.h>
#include <votca/ctp/qmpair.h>

namespace CSG = votca::csg;


namespace votca { namespace ctp {

class Topology;


class QMNBList : public CSG::PairList< Segment*, QMPair >
{
public:

    QMNBList(Topology* top) : _top(top), _cutoff(0) { };
   ~QMNBList() { CSG::PairList<Segment*, QMPair>::Cleanup(); }

    void    setCutoff(double cutoff) { _cutoff = cutoff; }
    double  getCutoff() { return _cutoff; }

    QMPair *Add(Segment* seg1, Segment* seg2);

    void PrintInfo(FILE *out);

protected:
    
    double      _cutoff;
    Topology   *_top;
};

}}


#endif	/* _QMNBList_H */

