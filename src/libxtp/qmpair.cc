/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#include <votca/xtp/qmpair.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/topology.h>

using namespace votca::tools;

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


  double QMPair::getLambdaO(int state)const {
    double result;
    if (state ==-1) result=_lambdaO_e;
    else if (state ==+1) result=_lambdaO_h;
    else if (state ==+2) result=_lambdaO_s;
    else if (state ==+3) result=_lambdaO_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
  }

  double QMPair::getRate12(int state) const{
    double result;
    if (state ==-1) result=_rate12_e;
    else if (state ==+1) result=_rate12_h;
    else if (state ==+2) result=_rate12_s;
    else if (state ==+3) result=_rate12_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
  }

  double QMPair::getRate21(int state) const{
    double result;
    if (state ==-1) result=_rate21_e;
    else if (state ==+1) result=_rate21_h;
    else if (state ==+2) result=_rate21_s;
    else if (state ==+3) result=_rate21_t;
    else throw std::runtime_error(" ERROR CODE whx__01l1o__");
    return result;
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

  bool QMPair::isPathCarrier(int carrier) const{
    bool result;
    if (carrier==-1) result=_has_e;
    else if (carrier==+1) result=_has_h;
    else if (carrier==+2) result=_has_s;
    else if (carrier==+3) result=_has_t;
    else throw std::runtime_error(" ERROR CODE whx__01p1r__");
    return result;
  }

  //only izindo uses this 
  void QMPair::setJs(const std::vector<double> Js, int state) {
    std::vector <double> ::const_iterator it;
    double Jeff2 = 0.0;
    if (state == -1) {
      for (it = Js.begin(); it < Js.end(); ++it) {
        Jeff2 += (*it)*(*it);
      }
      Jeff2 /= double(Js.size());
      _Jeff2_e = Jeff2; 
    }
    else if (state == +1) {
      for (it = Js.begin(); it < Js.end(); ++it) {
        Jeff2 += (*it)*(*it);
      }
      Jeff2 /= double(Js.size());
      _Jeff2_h = Jeff2;
    }
    else {
      throw std::runtime_error(" ERROR CODE whx__01x1u__");
    }
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

  double QMPair::getJeff2(int state) const{
    double result;
    if (state == -1) result =_Jeff2_e;
    else if (state == +1) result =_Jeff2_h;
    else if (state == +2) result =_Jeff2_s;
    else if (state == +3) result =_Jeff2_t;
    else throw std::runtime_error(" ERROR CODE whx__01s1j__");
    return result;
  }



  Segment *QMPair::Seg2PbCopy() {
    if (_hasGhost) {
      return _ghost;
    }
    else {
      return second;
    }
  }

}}
