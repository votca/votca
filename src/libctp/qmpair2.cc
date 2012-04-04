#include <votca/ctp/qmpair2.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

QMPair2::~QMPair2() {
    if (_ghost != NULL) {
        delete _ghost; _ghost = NULL;
    }
}


QMPair2::QMPair2(int id, Segment *seg1, Segment *seg2)
        : _id(id), std::pair<Segment*, Segment*>(seg1, seg2),
          _hasGhost(0),
          _rate12_e(0), _rate21_e(0),
          _rate12_h(0), _rate21_h(0),
          _has_e(false), _has_h(false),
          _lambdaO_e(0), _lambdaO_h(0),
          _Jeff2_e(0),   _Jeff2_h(0) {

    _top = seg1->getTopology();

    vec r1 = seg1->getPos();
    vec r2 = seg2->getPos();

    _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

    // Check whether pair formed across periodic boundary
    if ( abs(r2 - r1 - _R) > 1e-8 ) {

        _ghost = new Segment(seg2);
        _ghost->TranslateBy(r1 - r2 + _R);
        _hasGhost = true;

    }
    else {
        _ghost = NULL;
    }
}


double QMPair2::getLambdaO(int state) {

    return (state == -1) ? _lambdaO_e : _lambdaO_h;
}

double QMPair2::getRate12(int carrier) {

    return (carrier == -1) ? _rate12_e : _rate12_h;

}

double QMPair2::getRate21(int carrier) {

    return (carrier == -1) ? _rate21_e : _rate21_h;

}

void QMPair2::setLambdaO(double lO, int state) {

    if (state == -1) {
        _lambdaO_e = lO;
    }
    else if (state == +1) {
        _lambdaO_h = lO;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01l1o__");
    }
    

}

void QMPair2::setRate12(double rate, int state) {

    if (state == -1) {
        _rate12_e = rate;
    }
    else if (state == +1) {
        _rate12_h = rate;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01v1s__");
    }
}

void QMPair2::setRate21(double rate, int state) {

    if (state == -1) {
        _rate21_e = rate;
    }
    else if (state == +1) {
        _rate21_h = rate;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01w1t__");
    }
}

void QMPair2::setIsPathCarrier(bool yesno, int carrier) {

    if (carrier == -1) {
        _has_e = yesno;
    }
    else if (carrier == +1) {
        _has_h = yesno;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01p1r__");
    }
}

bool QMPair2::isPathCarrier(int carrier) {

    return (carrier == -1) ? _has_e : _has_h;
}


void QMPair2::setJs(const vector<double> Js, int state) {

    if (state == -1) {
        this->_Js_e = Js;
        this->calcJeff2(state);
    }
    else if (state == +1) {
        this->_Js_h = Js;
        this->calcJeff2(state);
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01x1u__");
    }
}


double QMPair2::calcJeff2(int state) {

    vector <double> ::iterator it;
    double Jeff2 = 0.0;

    if (state == -1) {
        for (it = _Js_e.begin(); it < _Js_e.end(); it++) {
            Jeff2 += (*it)*(*it);
        }
        Jeff2 /= double(_Js_e.size());
        _Jeff2_e = Jeff2;
    }

    else if (state == +1) {
        for (it = _Js_h.begin(); it < _Js_h.end(); it++) {
            Jeff2 += (*it)*(*it);
        }
        Jeff2 /= double(_Js_h.size());
        _Jeff2_h = Jeff2;
    }

    else {
        throw std::runtime_error(" ERROR CODE whx__01y1v__");
    }
    
    return Jeff2;
}

void QMPair2::setJeff2(double Jeff2, int state) {

    if (state == -1) {
        _Jeff2_e = Jeff2;
    }
    else if (state == +1) {
        _Jeff2_h = Jeff2;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01s1j__");
    }
}

Segment *QMPair2::Seg2PbCopy() {
    if (_hasGhost) {
        return _ghost;
    }
    else {
        return second;
    }
}

void QMPair2::WritePDB(string fileName) {

    FILE *pdb = NULL;
    pdb = fopen(fileName.c_str(), "w");

    this->first->WritePDB(pdb, "Atoms", "MD");
    this->second->WritePDB(pdb, "Atoms", "MD");
    if (this->HasGhost()) { this->_ghost->WritePDB(pdb, "Atoms", "MD"); }

    fclose(pdb);
}

}}
