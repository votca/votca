#ifndef __POLARSEG__H
#define	__POLARSEG__H

#include <votca/tools/vec.h>
#include <votca/ctp/apolarsite.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace votca { namespace ctp {

class PolarNb;
    
class PolarSeg : public vector<APolarSite*>
{

public:

    PolarSeg() : _id(-1), _is_charged(true), _is_polarizable(true) {};
    PolarSeg(int id, vector<APolarSite*> &psites);
    PolarSeg(PolarSeg *templ);
   ~PolarSeg();

    const int &getId() { return _id; }
    const vec &getPos() { return _pos; }
    
    void CalcPos();
    void CalcIsCharged();
    void CalcIsPolarizable();
    double CalcTotQ();    
    void Translate(const vec &shift);
    
    vector<PolarNb*> &PolarNbs() { return _nbs; }
    void ReservePolarNbs(int nbsize) { _nbs.reserve(nbsize); }
    void AddPolarNb(PolarNb *nb) { _nbs.push_back(nb); }
    void ClearPolarNbs();
    void PrintPolarNbPDB(string outfile);
    void WriteMPS(string mpsfile, string tag="");
    
    bool IsCharged() { return _is_charged; }    
    bool IsPolarizable() { return _is_polarizable; }
    
    template<class Archive>
    void serialize(Archive &arch, const unsigned int version) {
        arch & boost::serialization::base_object< vector<APolarSite*> >(*this);
        arch & _id;
        arch & _pos;
        arch & _is_charged;
        arch & _is_polarizable;
        return;
    }
    
private:

    int _id;
    vec _pos;
    bool _is_charged;
    bool _is_polarizable;
    vector<PolarNb*> _nbs;


};


class PolarNb
{
public:
    // _pbcshift = dr(nb,shift;pbc) - dr(nb,shift) + imageboxvector
    // Effective connection vector
    // top->ShortestConnect(ref,nb) + _pbcshift
    PolarNb(PolarSeg *nb, vec &dr12, vec &s22x) 
        : _nb(nb), _dr12(dr12), _s22x(s22x) {};        
    PolarSeg *getNb() { return _nb; }
    vec &getR() { return _dr12; }
    vec &getS() { return _s22x; }
private:
    PolarSeg *_nb;
    vec _dr12;
    vec _s22x;
};


}}

#endif