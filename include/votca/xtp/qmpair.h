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

#ifndef VOTCA_XTP_QMPAIR_H
#define VOTCA_XTP_QMPAIR_H

#include <vector>
#include <votca/tools/vec.h>
#include <utility>

#include <votca/xtp/segment.h>

namespace votca { namespace xtp {

class Topology;

class QMPair : public std::pair< Segment*, Segment* >
{
public:
    
    enum PairType{ 
        Hopping,
        Excitoncl,  
    };

    QMPair() :  _R(0,0,0),
                _ghost(NULL), 
                _top(NULL),
                _id(-1),   
                _hasGhost(false),
                _rate12_e(0),
                _rate21_e(0),
                _rate12_h(0),
                _rate21_h(0),
                _has_e(false),
                _has_h(false),
                _lambdaO_e(0),
                _lambdaO_h(0),
                _Jeff2_e(0),
                _Jeff2_h(0),
                _rate12_s(0),
                _rate21_s(0),
                _rate12_t(0),
                _rate21_t(0),
                _has_s(false),
                _has_t(false),
                _lambdaO_s(0),
                _lambdaO_t(0),   
                _Jeff2_s(0),
                _Jeff2_t(0),
                _pair_type(Hopping) { };
    QMPair(int id, Segment *seg1, Segment *seg2);
   ~QMPair();


   int       getId() { return _id; }
   void      setId(int id) { _id=id; }
   Topology *getTopology() { return _top; }
   void      setTopology(Topology *top) { _top = top; }
   votca::tools::vec      &R() { return _R; }
   double    Dist() { return abs(_R); }
   votca::tools::vec       getPos() { return 0.5*(first->getPos() + second->getPos()); }

   void     setIsPathCarrier(bool yesno, int carrier);
   bool     isPathCarrier(int carrier);

   void     setLambdaO(double lO, int carrier);
   double   getLambdaO(int carrier);
   
   double   getReorg12(int state) { return first->getU_nC_nN(state) + second->getU_cN_cC(state); } // 1->2
   double   getReorg21(int state) { return first->getU_cN_cC(state) + second->getU_nC_nN(state); } // 2->1
  
   double   getReorg12_x(int state) { return first->getU_nX_nN(state) + second->getU_xN_xX(state); } // 1->2
   double   getReorg21_x(int state) { return first->getU_xN_xX(state) + second->getU_nX_nN(state); } // 1->2

   void     setRate12(double rate, int state);
   void     setRate21(double rate, int state);
   double   getRate12(int state);
   double   getRate21(int state);
   votca::tools::vec      getR();

   //only used for compability reasons with izindo
   void     setJs(const std::vector <double> Js, int state);
   
   double   calcJeff2(int state);
   double   getJeff2(int state) ;
   void     setJeff2(double Jeff2, int state);
  

   double   getdE12(int state) { return second->getSiteEnergy(state)
                                       -first->getSiteEnergy(state); }

   Segment* Seg1PbCopy() { return first; }
   Segment* Seg2PbCopy();
   Segment* Seg1() { return first; }
   Segment* Seg2() { return second; }

   bool     HasGhost() { return _hasGhost; }

   // superexchange pairs have a list of bridging segments
   void     setType( PairType pair_type ) { _pair_type = pair_type; }
   void     setType( int pair_type ) { _pair_type = (PairType) pair_type; }
   PairType &getType(){return _pair_type;}

protected:

    votca::tools::vec         _R;

    Segment    *_ghost;
    Topology   *_top;
    int         _id;
    bool        _hasGhost;
    
    double _rate12_e;    // from ::Rates        output    DEFAULT 0
    double _rate21_e;    // from ::Rates        output    DEFAULT 0
    double _rate12_h;
    double _rate21_h;
    bool _has_e;       // from ::Rates        input     DEFAULT 0
    bool _has_h;
    double _lambdaO_e;   // from ::EOutersphere output    DEFAULT 0
    double _lambdaO_h;
    
    double          _Jeff2_e;
    double          _Jeff2_h;
    //excition part s:singlet t:triplet
    // state +2: singlet
    //state +3:triplet
    
    
    double _rate12_s;   
    double _rate21_s; 
    double _rate12_t;
    double _rate21_t;
    bool _has_s;       
    bool _has_t;
    double _lambdaO_s;   
    double _lambdaO_t; 
    double          _Jeff2_s;
    double          _Jeff2_t;

    PairType _pair_type;
    
};

}}


#endif // VOTCA_XTP_QMPAIR_H
