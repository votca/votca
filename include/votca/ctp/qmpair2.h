#ifndef _QMPAIR2_H
#define _QMPAIR2_H

#include "segment.h"
#include <utility>


namespace votca { namespace ctp {

class Topology;



class QMPair2 : public std::pair< Segment*, Segment* >
{
public:
    QMPair2() : _R(0,0,0), _rate_12(0), _rate_21(0), _ghost(NULL),
                _top(NULL), _id(-1), _hasGhost(0) { };
    QMPair2(int id, Segment *seg1, Segment *seg2);
   ~QMPair2() { if (_ghost != NULL) { delete _ghost; _ghost = NULL; } };


   int      getId() { return _id; }
   vec     &R() { return _R; }
   double   Dist() { return abs(_R); }

   void     setRate12(double rate) { _rate_12 = rate; }
   void     setRate21(double rate) { _rate_21 = rate; }

   void     setJs(const vector <double> &Js) { _Js = Js; }
   double   calcJeff2();
   vector<double> &Js() { return _Js; }

   Segment* Seg1PbCopy() { return first; }
   Segment* Seg2PbCopy();

   bool     HasGhost() { return _hasGhost; }



protected:

    int         _id;
    vec         _R;
    Topology   *_top;

    Segment    *_ghost;
    bool        _hasGhost;

    double      _rate_12;
    double      _rate_21;
    vector <double> _Js;




};

}}


#endif
