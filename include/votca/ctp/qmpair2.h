#ifndef _QMPAIR2_H
#define _QMPAIR2_H

#include "segment.h"
#include <utility>


namespace votca { namespace ctp {

class Topology;



class QMPair2 : public std::pair< Segment*, Segment* >
{
public:
    QMPair2() : _R(0,0,0), _ghost(NULL), _top(NULL),
                _id(-1), _hasGhost(0) { };
    QMPair2(int id, Segment *seg1, Segment *seg2);
   ~QMPair2();


   int       getId() { return _id; }
   Topology *getTopology() { return _top; }
   void      setTopology(Topology *top) { _top = top; }
   vec      &R() { return _R; }
   double    Dist() { return abs(_R); }

   void     setLambdaO(int carrier, double lbd);
   void     setLambdaO(double lambda) { _lambdaO = lambda; }
   bool     hasLambdaO() { return _hasLambdaO; }
   double  &getLambdaO(int carrier);
   double  &getLambdaO() { return _lambdaO; }

   void     setRate12(int carrier, double rate);
   void     setRate12(double rate) { _rate12 = rate; }
   void     setRate21(int carrier, double rate);
   void     setRate21(double rate) { _rate21 = rate; }
   bool     hasRates() { return _hasRates; }
   double  &getRate12(int carrier);
   double  &getRate12() { return _rate12; }
   double  &getRate21(int carrier);
   double  &getRate21() { return _rate21; }

   void     setJs(const vector <double> Js) { _Js = Js; }
   double   calcJeff2();
   vector<double> &Js() { return _Js; }

   Segment* Seg1PbCopy() { return first; }
   Segment* Seg2PbCopy();
   Segment* Seg1() { return first; }
   Segment* Seg2() { return second; }

   bool     HasGhost() { return _hasGhost; }



protected:

    int         _id;
    vec         _R;
    Topology   *_top;

    Segment    *_ghost;
    bool        _hasGhost;

    double _lambdaO;
    double _rate12;
    double _rate21;

    map< int,       double > _lambdasO;
    //   -1(=> e)   lambda for electrons
    //   +1(=> h)   lambda for holes
    bool _hasLambdaO;

    map< int,       double > _rates12;
    map< int,       double > _rates21;
    //   -1(=> e)   rate for electrons
    //   +1(=> h)   rate for holes
    bool _hasRates;

    map< int,       bool > _hasCarrier; // not used right now
    //   -1(=> e)   true or false
    //   +1(=> h)   true or false
    
    vector <double> _Js;




};

}}


#endif
