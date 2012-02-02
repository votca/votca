#ifndef _QMNBLIST2_H
#define	_QMNBLIST2_H



#include <votca/csg/pairlist.h>
#include <votca/ctp/qmpair2.h>

namespace CSG = votca::csg;


namespace votca { namespace ctp {

class Topology;


class QMNBList2 : public CSG::PairList< Segment*, QMPair2 >
{
public:

    QMNBList2(Topology* top) : _top(top), _cutoff(0) { };
   ~QMNBList2() { CSG::PairList<Segment*, QMPair2>::Cleanup(); }

    void    setCutoff(double cutoff) { _cutoff = cutoff; }
    double  getCutoff() { return _cutoff; }

    QMPair2 *Add(Segment* seg1, Segment* seg2);

    void PrintInfo(FILE *out);

protected:
    
    double      _cutoff;
    Topology   *_top;
};

}}


#endif	/* _QMNBLIST2_H */

