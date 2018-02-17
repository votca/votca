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


#ifndef _VOTCA_XTP_POLARSITE_H
#define _VOTCA_XTP_POLARSITE_H

#include <vector>
#include <string>

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>

namespace votca { namespace xtp {

class Topology;
class Molecule;
class Segment;
class Fragment;


class PolarSite
{

    friend class EMultipole;
    friend class EMultipole_StdAl;
    friend class EOutersphere;
    friend class ECoulomb;
    friend class Interactor;
    friend class Interactor3;
    friend class InteractorMod;
    friend class XMP;
    friend class XInteractor;

public:

    PolarSite(int id, std::string name)
            : _id(id), _name(name), _locX(votca::tools::vec(1,0,0)),
              _locY(votca::tools::vec(0,1,0)),    _locZ(votca::tools::vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

    PolarSite() 
            : _id(-1),  _locX(votca::tools::vec(1,0,0)),
              _locY(votca::tools::vec(0,1,0)), _locZ(votca::tools::vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

   ~PolarSite() {};

    int             &getId() { return _id; }
    std::string          &getName() { return _name; }
    votca::tools::vec             &getPos() { return _pos; }
    int             &getRank() { return _rank; }
    Topology        *getTopology() { return _top; }
    Segment         *getSegment() { return _seg; }
    Fragment        *getFragment() { return _frag; }

    void            setPos(votca::tools::vec &pos) { _pos = pos; }
    void            setRank(int rank) { _rank = rank; } // rank; } // OVERRIDE
    void            setTopology(Topology *top) { _top = top; }
    void            setSegment(Segment *seg) { _seg = seg; }
    void            setFragment(Fragment *frag) { _frag = frag; }

    std::vector<double> &getQs(int state) { return _Qs[state+1]; }
    void            setQs(std::vector<double> Qs, int state) { _Qs[state+1] = Qs; }
    void            setPs(double polar, int state) { _Ps[state+1] = polar; }
    double         &getPs(int state) { return _Ps[state+1]; }
    double         &getP1() { return P1; }
    double         &getQ00() { return Q00; }
    void            Charge(int state);
    void            ChargeDelta(int state1, int state2);

    void            Induce(double wSOR = 0.25);
    void            InduceDirect();
    void            ResetFieldU() { FUx = FUy = FUz = 0.0; }
    void            ResetFieldP() { FPx = FPy = FPz = 0.0; }
    void            ResetU1Hist() { U1_Hist.clear(); }
    void            Depolarize();
    double          HistdU();


    void            ImportFrom(PolarSite *templ, std::string tag = "basic");
    void            Translate(const votca::tools::vec &shift);
    void            Rotate(const votca::tools::matrix &rot, const votca::tools::vec &refPos);

    void            PrintInfo(std::ostream &out);
    void            PrintInfoInduce(std::ostream &out);
    void            PrintInfoVisual(FILE *out);
    void            PrintPDB(FILE *out, votca::tools::vec shift);
    void            WriteChkLine(FILE *, votca::tools::vec &, bool, std::string, double);
    void            WriteXyzLine(FILE *, votca::tools::vec &, std::string);




private:

    int     _id;
    std::string  _name;
    votca::tools::vec     _pos;
    votca::tools::vec     _locX;
    votca::tools::vec     _locY;
    votca::tools::vec     _locZ;

    Topology *_top;
    Segment  *_seg;
    Fragment *_frag;
    
    std::vector < std::vector<double> > _Qs;
    int     _rank;
                           
    std::vector < double > _Ps;
    double P1;                              // Dipole polarizability

    double Q00;
    double Q1x, Q1y, Q1z;    
    double Q20, Q21c, Q21s, Q22c, Q22s;

    double U1x, U1y, U1z;                   // Induced dipole
    double FPx, FPy, FPz;                   // Electric field (due to permanent)
    double FUx, FUy, FUz;                   // Electric field (due to induced)
    std::vector< votca::tools::vec > U1_Hist;                  // Ind. u history



    

};


}}

#endif
