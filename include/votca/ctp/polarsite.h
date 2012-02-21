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

public:

    PolarSite(int id, string name) : _id(id), _name(name) { _Qs.resize(3); };
   ~PolarSite() {};


    vec &getPos() { return _pos; }
    int &getRank() { return _rank; }
    

    void setPos(vec &pos) { _pos = pos; }
    void setRank(int rank) { _rank = rank; }

    // Multipole management of charge states
    vector<double> &getQs(int state) { return _Qs[state+1]; }
    void            setQs(vector<double> Qs, int state) { _Qs[state+1] = Qs; }

    void setTopology(Topology *top) { _top = top; }
    void setSegment(Segment *seg) { _seg = seg; }
    void setFragment(Fragment *frag) { _frag = frag; }

    void ImportFrom(PolarSite *templ);


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

    int     _rank;
    vector < vector<double> > _Qs;

    

};






}}

#endif