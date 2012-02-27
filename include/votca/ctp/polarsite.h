#ifndef POLARSITE_H
#define POLARSITE_H

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>

namespace votca { namespace ctp {

using namespace votca::tools;


class Topology;
class Molecule;
class Segment;
class Fragment;


class PolarSite
{

    friend class EMultipole2;
    friend class Interactor;

public:

    PolarSite(int id, string name)
            : _id(id), _name(name), _locX(vec(1,0,0)),
              _locY(vec(0,1,0)),    _locZ(vec(0,0,1)) { _Qs.resize(3); };

    PolarSite() 
            : _id(-1),  _locX(vec(1,0,0)),
              _locY(vec(0,1,0)), _locZ(vec(0,0,1))    { _Qs.resize(3); };

   ~PolarSite() {};

    int             &getId() { return _id; }
    string          &getName() { return _name; }
    vec             &getPos() { return _pos; }
    int             &getRank() { return _rank; }
    Topology        *getTopology() { return _top; }
    Segment         *getSegment() { return _seg; }
    Fragment        *getFragment() { return _frag; }

    void            setPos(vec &pos) { _pos = pos; }
    void            setRank(int rank) { _rank = 1; } // rank; } // OVERRIDE
    void            setTopology(Topology *top) { _top = top; }
    void            setSegment(Segment *seg) { _seg = seg; }
    void            setFragment(Fragment *frag) { _frag = frag; }

    vector<double> &getQs(int state) { return _Qs[state+1]; }
    void            setQs(vector<double> Qs, int state) { _Qs[state+1] = Qs; }
    void            Charge(int state);
    void            setAlpha(double polarity) { alpha = polarity; }
    double         &getAlpha() { return alpha; }
    void            Depolarize() { U1x = U1y = U1z = 0.0; U1_Hist.clear(); }
    void            Add2IndHist(vec &u) { U1_Hist.push_back(u); }


    void            ImportFrom(PolarSite *templ, string tag = "basic");
    void            Translate(const vec &shift);
    void            Rotate(const matrix &rot, const vec &refPos);

    void            PrintInfo(std::ostream &out);




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

    double alpha;                           // Polarizability

    double Q00;
    double Q1x, Q1y, Q1z;
    double U1x, U1y, U1z;                   // Induced dipole
    double Q20, Q21c, Q21s, Q22c, Q22s;

    double Fx, Fy, Fz;                      // Electric field
    vector< vec > U1_Hist;                  // Ind. u history



    

};






}}

#endif