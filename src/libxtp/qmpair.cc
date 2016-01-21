/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#include <votca/xtp/qmpair.h>
#include <votca/xtp/topology.h>

namespace votca { namespace xtp {

QMPair::~QMPair() {
    if (_ghost != NULL) {
        delete _ghost; _ghost = NULL;
    }
}


QMPair::QMPair(int id, Segment *seg1, Segment *seg2)
        : std::pair<Segment*, Segment*>(seg1, seg2), _id(id),
          _hasGhost(0),
          _rate12_e(0), _rate21_e(0),
          _rate12_h(0), _rate21_h(0),      
          _has_e(false), _has_h(false),
          _lambdaO_e(0), _lambdaO_h(0),
         _Jeff2_e(0),   _Jeff2_h(0),
          _rate12_s(0), _rate21_s(0),
          _rate12_t(0), _rate21_t(0),
          _has_s(false), _has_t(false),       
          _lambdaO_s(0), _lambdaO_t(0),
          _Jeff2_s(0),   _Jeff2_t(0),_pair_type( Hopping ) {

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
    double result;
    if (state ==-1) result=_lambdaO_e;
    else if (state ==+1) result=_lambdaO_h;
    else if (state ==+2) result=_lambdaO_s;
    else if (state ==+3) result=_lambdaO_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
}

double QMPair::getRate12(int state) {
    double result;
    if (state ==-1) result=_rate12_e;
    else if (state ==+1) result=_rate12_h;
    else if (state ==+2) result=_rate12_s;
    else if (state ==+3) result=_rate12_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
}

double QMPair::getRate21(int state) {
    double result;
    if (state ==-1) result=_rate21_e;
    else if (state ==+1) result=_rate21_h;
    else if (state ==+2) result=_rate21_s;
    else if (state ==+3) result=_rate21_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
}

vec QMPair::getR() {

    return _R;
}


void QMPair::setLambdaO(double lO, int state) {  
    if (state ==-1) _lambdaO_e = lO;
    else if (state ==+1) _lambdaO_h = lO;
    else if (state ==+2) _lambdaO_s = lO;
    else if (state ==+3) _lambdaO_t = lO;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
}

void QMPair::setRate12(double rate, int state) {
    if (state ==-1) _rate12_e = rate;
    else if (state ==+1) _rate12_h = rate;
    else if (state ==+2) _rate12_s = rate;
    else if (state ==+3) _rate12_t = rate;
    else throw std::runtime_error(" ERROR CODE whx__01v1s__");
}

void QMPair::setRate21(double rate, int state) {
    if (state ==-1) _rate21_e = rate;
    else if (state ==+1) _rate21_h = rate;
    else if (state ==+2) _rate21_s = rate;
    else if (state ==+3) _rate21_t = rate;
    else throw std::runtime_error(" ERROR CODE whx__01w1t__");
}

void QMPair::setIsPathCarrier(bool yesno, int carrier) {
    if (carrier == -1)_has_e = yesno;
    else if (carrier == +1)_has_h = yesno;
    else if (carrier == +2)_has_s = yesno;
    else if (carrier == +3)_has_t = yesno;
    else throw std::runtime_error(" ERROR CODE whx__01p1r__");
}

bool QMPair::isPathCarrier(int carrier) {
    bool result;
    if (carrier==-1) result=_has_e;
    else if (carrier==+1) result=_has_h;
    else if (carrier==+2) result=_has_s;
    else if (carrier==+3) result=_has_t;
    else throw std::runtime_error(" ERROR CODE whx__01p1r__");
    return result;
}


void QMPair::setJs(const vector<double> Js, int state) {

    if (state == -1) {
        this->_Js_e = Js;
        //double Jeff2 =
	(void)this->calcJeff2(state);
    }
    else if (state == +1) {
        this->_Js_h = Js;
        //double Jeff2 =
	(void)this->calcJeff2(state);
    }
    else if (state == +2) {
        this->_Js_s = Js;
        (void)this->calcJeff2(state);
    }
    else if (state == +3) {
        this->_Js_t = Js;
        (void)this->calcJeff2(state);
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
    
    else if (state == +2) {
        for (it = _Js_s.begin(); it < _Js_s.end(); it++) {
            Jeff2 += (*it)*(*it);
        }
        Jeff2 /= double(_Js_s.size());
        _Jeff2_s = Jeff2;
    }
    
    else if (state == +3) {
        for (it = _Js_t.begin(); it < _Js_t.end(); it++) {
            Jeff2 += (*it)*(*it);
        }
        Jeff2 /= double(_Js_t.size());
        _Jeff2_t = Jeff2;
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
     else if (state == +2) {
        _Jeff2_s = Jeff2;
    }
     else if (state == +3) {
        _Jeff2_t = Jeff2;
    }
    else {
        throw std::runtime_error(" ERROR CODE whx__01s1j__");
    }
}

 double QMPair::getJeff2(int state) {
     double result;
     if (state == -1) result =_Jeff2_e;
     else if (state == +1) result =_Jeff2_h;
     else if (state == +2) result =_Jeff2_s;
     else if (state == +3) result =_Jeff2_t;
     else throw std::runtime_error(" ERROR CODE whx__01s1j__");
     return result;
}
 
 vector<double> &QMPair::Js(int state) {
     if (state==-1) return _Js_e;
     else if (state==+1) return _Js_h;
     else if (state==+2) return _Js_s;
     else if (state==+3) return _Js_t;
     else throw std::runtime_error(" ERROR CODE whx__01s1j__");
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
