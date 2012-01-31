#include <votca/ctp/qmpair2.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

QMPair2::QMPair2(int id, Segment *seg1, Segment *seg2)
        : _id(id), std::pair<Segment*, Segment*>(seg1, seg2),
          _rate_12(0), _rate_21(0), _hasGhost(0) {

    _top = seg1->getTopology();

    vec r1 = seg1->getPos();
    vec r2 = seg2->getPos();

    _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

    // Check whether pair formed across periodic boundary
    if ( abs(r2 - r1 - _R) > 1e-8 ) {
        _ghost = new Segment( second->getId(), second->getName() );
        _ghost->setPos(r1 + _R);
        _ghost->setMolecule( second->getMolecule() );
        _ghost->setTopology( _top );
        _hasGhost = 1;
    }
    else {
        _ghost = NULL;
    }
}


double QMPair2::calcJeff2() {
    vector <double> ::iterator it;
    double Jeff2 = 0;
    for (it = _Js.begin(); it < _Js.end(); it++) {
        Jeff2 += (*it)*(*it);
    }
    Jeff2 /= double(_Js.size());
    return Jeff2;
}

Segment *QMPair2::Seg2PbCopy() {
    if (_hasGhost) {
        return _ghost;
    }
    else {
        return second;
    }
}

}}
