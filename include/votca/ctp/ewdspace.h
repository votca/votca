#ifndef VOTCA_CTP_EWDSPACE_H
#define VOTCA_CTP_EWDSPACE_H

#include <cmath>
#include <votca/tools/vec.h>

namespace votca {
namespace ctp {
namespace EWD {

// DEFINES THE FOLLOWING
// o *cmplx* structure
// o *triple* (templated on contained type)
// o *VectorSort* (functor, templated on type and norm)
//    - *MaxNorm*
//    - *EucNorm*
//    - *KNorm*
// o KVector (with 'grade' attribute)
// o Conversion constants to harmonize with EwdInteractor
    
    
// "static" not really required here
static const double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
static const double int2V_m = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-18;
static const double rSqrtPi = 0.564189583547756279280349644977832;

struct cmplx
{
    cmplx() { ; }
    cmplx(double re, double im) : _re(re), _im(im) { ; }
    cmplx(const cmplx &c) : _re(c._re), _im(c._im) { ; }
    cmplx &operator*=(const double &d);
    cmplx &operator+=(const cmplx &c);
    cmplx &operator-=(const cmplx &c);
    cmplx &Conjugate() { _im *= -1; return *this; }
    double _re;
    double _im;
};

inline cmplx &cmplx::operator*=(const double &d) { 
    _re*=d;
    _im*=d;
    return *this;
}

inline cmplx &cmplx::operator+=(const cmplx &c) {
    _re+=c._re;
    _im+=c._im;
    return *this;
}

inline cmplx &cmplx::operator-=(const cmplx &c) {
    _re-=c._re;
    _im-=c._im;
    return *this;
}

inline cmplx operator+(const cmplx &cl, const cmplx &cr) {
    return (cmplx(cl)+=cr);
}

inline cmplx operator-(const cmplx &cl, const cmplx &cr) {
    return (cmplx(cl)-=cr);
}

inline cmplx operator*(const cmplx &cl, const double &d) {
    return (cmplx(cl)*=d);
}

inline cmplx operator*(const double &d, const cmplx &cr) {
    return (cmplx(cr)*=d);
}

    
template<typename NrTyp=double>
struct triple
{
    triple() { ; }
    triple(NrTyp pp, NrTyp pu, NrTyp uu) : _pp(pp), _pu(pu), _uu(uu) { ; }
    triple(const triple<NrTyp> &t) : _pp(t._pp), _pu(t._pu), _uu(t._uu) { ; }
    triple &operator*=(const double &d);
    triple &operator+=(const triple &t);
    triple &operator-=(const triple &t);
    NrTyp Sum() { return _pp + _pu + _uu; }
    NrTyp _pp;
    NrTyp _pu;
    NrTyp _uu;
};

template<typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator*=(const double &d) {
    _pp*=d;
    _pu*=d;
    _uu*=d;
    return *this;
}

template<typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator+=(const triple<NrTyp> &t) {
    _pp+=t._pp;
    _pu+=t._pu;
    _uu+=t._uu;
    return *this;
}

template<typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator-=(const triple<NrTyp> &t) {
    _pp-=t._pp;
    _pu-=t._pu;
    _uu-=t._uu;
    return *this;
}

template<typename NrTyp>
inline triple<NrTyp> operator+(const triple<NrTyp> &tl, const triple<NrTyp> &tr) {
    return (triple<NrTyp>(tl)+=tr);
}

template<typename NrTyp>
inline triple<NrTyp> operator-(const triple<NrTyp> &tl, const triple<NrTyp> &tr) {
    return (triple<NrTyp>(tl)-=tr);
}

template<typename NrTyp>
inline triple<NrTyp> operator*(const triple<NrTyp> &tl, const double &d) {
    return (triple<NrTyp>(tl)*=d);
}

template<typename NrTyp>
inline triple<NrTyp> operator*(const double &d, const triple<NrTyp> &tr) {
    return (triple<NrTyp>(tr)*=d);
}


// To sort K-vectors via std::sort using a norm functor
template<class Norm, class V>
struct VectorSort
{
    VectorSort() : _p(1e-40) { ; }
    VectorSort(double precision) : _p(precision) { ; }
    inline bool operator() (const V &v1, const V &v2);
    inline bool MatchDouble(double a, double b) 
        { return ((a-b)*(a-b) < _p) ? true : false; }
    double _p;
    Norm _norm;
};

// Tschebyschow norm functor
struct MaxNorm { inline double operator() (const votca::tools::vec &v) 
    { return votca::tools::maxnorm(v); } };
// Euclidean norm functor
struct EucNorm { inline double operator() (const votca::tools::vec &v) 
    { return votca::tools::abs(v); } };

// K-vector class (for grading purposes)
struct KVector
{
    KVector(votca::tools::vec k, double grade)
        : _k(k), _grade(grade) { ; }

    votca::tools::vec _k;
    double _grade;

    const votca::tools::vec &getK() const { return _k; }
    const double &getGrade() const { return _grade; }
    const double &getX() const { return _k.getX(); }
    const double &getY() const { return _k.getY(); }
    const double &getZ() const { return _k.getZ(); }            
};

// Specialized K-vector norm
struct KNorm { inline double operator() (const KVector &v)
    { return -v.getGrade(); } };

    
template<class Norm, class V>
inline bool VectorSort<Norm,V>::operator() (const V &v1,
    const V &v2) {
    bool smaller = false;
    // LEVEL 1: MAGNITUDE
    double V1 = _norm(v1);
    double V2 = _norm(v2);
    if (MatchDouble(V1,V2)) {
        // LEVEL 2: X
        double X1 = v1.getX();
        double X2 = v2.getX();
        if (MatchDouble(X1,X2)) {
            // LEVEL 3: Y
            double Y1 = v1.getY();
            double Y2 = v2.getY();
            if (MatchDouble(Y1,Y2)) {
                // LEVEL 4: Z
                double Z1 = v1.getZ();
                double Z2 = v2.getZ();
                if (MatchDouble(Z1,Z2)) smaller = true;
                else smaller = (Z1 < Z2) ? true : false;
            }
            else smaller = (Y1 < Y2) ? true : false;
        }
        else smaller = (X1 < X2) ? true : false;
    }
    else smaller = (V1 < V2) ? true : false;          
    return smaller;
}


}}}

#endif