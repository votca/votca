#include <votca/ctp/qmpair2.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

QMPair2::~QMPair2() {
    if (_ghost != NULL) {
        delete _ghost; _ghost = NULL;
    }
    _lambdasO.clear();
    _rates12.clear();
    _rates21.clear();
    _hasCarrier.clear();
}


QMPair2::QMPair2(int id, Segment *seg1, Segment *seg2)
        : _id(id), std::pair<Segment*, Segment*>(seg1, seg2),
          _hasGhost(0), _lambdaO(0), _rate12(0), _rate21(0) {

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


double &QMPair2::getLambdaO(int carrier) {
    return _lambdasO.at(carrier);
}

double &QMPair2::getRate12(int carrier) {
    return _rates12.at(carrier);
}

double &QMPair2::getRate21(int carrier) {
    return _rates21.at(carrier);
}

void QMPair2::setLambdaO(int carrier, double lbd) {
    _lambdasO[carrier] = lbd;
    _hasLambdaO = true;
}

void QMPair2::setRate12(int carrier, double rate) {
    _rates12[carrier] = rate;
    _hasRates = true;
}

void QMPair2::setRate21(int carrier, double rate) {
    _rates21[carrier] = rate;
    _hasRates = true;
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
