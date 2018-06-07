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

#ifndef __VOTCA_XTP_POLARSITE_H
#define __VOTCA_XTP_POLARSITE_H


#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>

namespace votca { namespace xtp {
    
class PolarSite
{

public:

    PolarSite(int id, std::string name)
            : _id(id),              _name(name),         _isVirtual(false), 
              _locX(votca::tools::vec(1,0,0)),    _locY(votca::tools::vec(0,1,0)),   _locZ(votca::tools::vec(0,0,1)), 
              _top(0),              _seg(0),             _frag(0),
              _resolution(atomistic),PhiP(0.0),          PhiU(0.0)
            { _Qs.resize(3); _Ps.resize(3); this->Depolarize();
              for (int s = -1; s < 2; ++s) _Ps[s+1].ZeroMatrix(); }
    APolarSite(APolarSite *templ);
   ~APolarSite() {};
   
  
    
    // GET & SET & IMPORT FUNCTIONS
    int                 &getId() { return _id; }
    std::string         &getName() { return _name; }
    votca::tools::vec   &getPos() { return _pos; }
    // Point charge 0,1 dipole,2 quad 
    int            &getRank() { return _rank; }
    Topology       *getTopology() { return _top; }
    Segment        *getSegment() { return _seg; }
    Fragment       *getFragment() { return _frag; }
    // Don't know

    void            setTopology(Topology *top) { _top = top; }
    void            setSegment(Segment *seg) { _seg = seg; }
    void            setFragment(Fragment *frag) { _frag = frag; }    
    
    // COORDINATE TRANSFORMATION
    void            Translate(const votca::tools::vec &shift);
    void            Rotate(const votca::tools::matrix &rot, const votca::tools::vec &refPos);
 
    Eigen::VectorXd getQuad_cartesian();
    Eigen::VectorXd getQuad();
    Eigen::VectorXd getDipole();
    
    const Eigen::VectorXd & getMultipoles()const {return _multipoles;}


private:

    int     _id;
    std::string  _name;
    tools::vec     _pos;
    
    Topology *_top;
    Segment  *_seg;
    Fragment *_frag;
    
    int     _rank;

    Eigen::Matrix3d _Ps;
    Eigen::Vector3d _localpermanetField;
    Eigen::Vector3d _localinducedField;
    Eigen::Vector3d _inducedDipole;
    Eigen::Vector3d _inducedDipole_old;
    
    Eigen::VectorXd _multipoles; //Q00,Q1x, Q1y, Q1z,Q20, Q21c, Q21s, Q22c, Q22s,
    
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)
    
    
    // These values are basically the same as what is in the _Qs but in a different
    // format apparently for performance reasons... 
    double Q00;
    double Q1x, Q1y, Q1z;
    double Q20, Q21c, Q21s, Q22c, Q22s;
    double Qxx, Qxy, Qxz, Qyy, Qyz, Qzz;

    double U1x, U1y, U1z;                   // Induced dipole
    double FPx, FPy, FPz;                   // Electric field (due to permanent)
    double FUx, FUy, FUz;                   // Electric field (due to induced)
    std::vector< votca::tools::vec > U1_Hist;                  // Ind. u history
    res_t   _resolution;
   
   


};


}}


#endif
