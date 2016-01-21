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


#ifndef POLARSITE_H
#define POLARSITE_H

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>

namespace votca { namespace xtp {

using namespace votca::tools;


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

    PolarSite(int id, string name)
            : _id(id), _name(name), _locX(vec(1,0,0)),
              _locY(vec(0,1,0)),    _locZ(vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

    PolarSite() 
            : _id(-1),  _locX(vec(1,0,0)),
              _locY(vec(0,1,0)), _locZ(vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

   ~PolarSite() {};

    int             &getId() { return _id; }
    string          &getName() { return _name; }
    vec             &getPos() { return _pos; }
    int             &getRank() { return _rank; }
    Topology        *getTopology() { return _top; }
    Segment         *getSegment() { return _seg; }
    Fragment        *getFragment() { return _frag; }

    void            setPos(vec &pos) { _pos = pos; }
    void            setRank(int rank) { _rank = rank; } // rank; } // OVERRIDE
    void            setTopology(Topology *top) { _top = top; }
    void            setSegment(Segment *seg) { _seg = seg; }
    void            setFragment(Fragment *frag) { _frag = frag; }

    vector<double> &getQs(int state) { return _Qs[state+1]; }
    void            setQs(vector<double> Qs, int state) { _Qs[state+1] = Qs; }
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


    void            ImportFrom(PolarSite *templ, string tag = "basic");
    void            Translate(const vec &shift);
    void            Rotate(const matrix &rot, const vec &refPos);

    void            PrintInfo(std::ostream &out);
    void            PrintInfoInduce(std::ostream &out);
    void            PrintInfoVisual(FILE *out);
    void            PrintPDB(FILE *out, vec shift);
    void            WriteChkLine(FILE *, vec &, bool, string, double);
    void            WriteXyzLine(FILE *, vec &, string);




private:

    int     _id;
    string  _name;
    vec     _pos;
    vec     _locX;
    vec     _locY;
    vec     _locZ;

    Topology *_top;
    Segment  *_seg;
    Fragment *_frag;
    
    vector < vector<double> > _Qs;
    int     _rank;
                           
    vector < double > _Ps;
    double P1;                              // Dipole polarizability

    double Q00;
    double Q1x, Q1y, Q1z;    
    double Q20, Q21c, Q21s, Q22c, Q22s;

    double U1x, U1y, U1z;                   // Induced dipole
    double FPx, FPy, FPz;                   // Electric field (due to permanent)
    double FUx, FUy, FUz;                   // Electric field (due to induced)
    vector< vec > U1_Hist;                  // Ind. u history



    

};






}}

#endif
