/*
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/ctp/qmpair.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {

QMPair::~QMPair() {
    if (_ghost != NULL) {
        delete _ghost; _ghost = NULL;
    }
}


QMPair::QMPair(int id, Segment *seg1, Segment *seg2)
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


double QMPair::getLambdaO(int state) {

    return (state == -1) ? _lambdaO_e : _lambdaO_h;
}

double QMPair::getRate12(int carrier) {

    return (carrier == -1) ? _rate12_e : _rate12_h;

}

double QMPair::getRate21(int carrier) {

    return (carrier == -1) ? _rate21_e : _rate21_h;

}

void QMPair::setLambdaO(double lO, int state) {

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

void QMPair::setRate12(double rate, int state) {

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

void QMPair::setRate21(double rate, int state) {

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

void QMPair::setIsPathCarrier(bool yesno, int carrier) {

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

bool QMPair::isPathCarrier(int carrier) {

    return (carrier == -1) ? _has_e : _has_h;
}


void QMPair::setJs(const vector<double> Js, int state) {

    if (state == -1) {
        this->_Js_e = Js;
        double Jeff2 = this->calcJeff2(state);
    }
    else if (state == +1) {
        this->_Js_h = Js;
        double Jeff2 = this->calcJeff2(state);
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01x1u__");
    }
}


double QMPair::calcJeff2(int state) {

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

void QMPair::setJeff2(double Jeff2, int state) {

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

Segment *QMPair::Seg2PbCopy() {
    if (_hasGhost) {
        return _ghost;
    }
    else {
        return second;
    }
}

void QMPair::WritePDB(string fileName) {

    FILE *pdb = NULL;
    pdb = fopen(fileName.c_str(), "w");

    this->first->WritePDB(pdb, "Atoms", "MD");
    this->second->WritePDB(pdb, "Atoms", "MD");
    if (this->HasGhost()) { this->_ghost->WritePDB(pdb, "Atoms", "MD"); }

    fclose(pdb);
}

void QMPair::WriteXYZ(FILE *out, bool useQMPos) {

    int qmatoms = 0;

    vector< Atom* > ::iterator ait;

    for (ait = Seg1PbCopy()->Atoms().begin();
         ait < Seg1PbCopy()->Atoms().end();
         ++ait) {

        if ((*ait)->HasQMPart() || !useQMPos) {
            ++qmatoms;
        }
    }

    for (ait = Seg2PbCopy()->Atoms().begin();
         ait < Seg2PbCopy()->Atoms().end();
         ++ait) {

        if ((*ait)->HasQMPart() || !useQMPos) {
            ++qmatoms;
        }
    }

    fprintf(out, "%6d \n", qmatoms);
    fprintf(out, "\n");

    for (ait = Seg1PbCopy()->Atoms().begin();
         ait < Seg1PbCopy()->Atoms().end();
         ++ait) {

        if (!(*ait)->HasQMPart() && useQMPos) {
            continue;
        }

        vec pos;
        if (useQMPos) pos = (*ait)->getQMPos();
        else pos = (*ait)->getPos();
        
        string  name = (*ait)->getElement();

        fprintf(out, "%2s %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10);
    }

    for (ait = Seg2PbCopy()->Atoms().begin();
         ait < Seg2PbCopy()->Atoms().end();
         ++ait) {

        if (!(*ait)->HasQMPart() && useQMPos) {
            continue;
        }

        vec pos;
        if (useQMPos) pos = (*ait)->getQMPos();
        else pos = (*ait)->getPos();
        
        string  name = (*ait)->getElement();

        fprintf(out, "%2s %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10);
    }
}


}}
