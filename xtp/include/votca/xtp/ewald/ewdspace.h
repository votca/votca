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
#ifndef VOTCA_XTP_EWDSPACE_H
#define VOTCA_XTP_EWDSPACE_H

#include <cmath>
// #include <votca/tools/vec.h>
#include <boost/format.hpp>
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {
typedef Eigen::Vector3d vec;
typedef Eigen::Matrix3d matrix;

namespace EWD {

using namespace votca::tools::conv;
// DEFINES THE FOLLOWING
// o *cmplx* structure
// o *triple* (templated on contained type)
// o *VectorSort* (functor, templated on type and norm)
//    - *MaxNorm*
//    - *EucNorm*
//    - *KNorm*
// o KVector (with 'grade' attribute)
// o Conversion constants to blend well with EwdInteractor

struct cmplx {
  cmplx() { ; }
  cmplx(double re, double im) : _re(re), _im(im) { ; }
  // cmplx(const cmplx &c) : _re(c._re), _im(c._im) { ; }
  const double &Re() const { return _re; }
  const double &Im() const { return _im; }
  cmplx &operator*=(const double &d);
  cmplx &operator+=(const cmplx &c);
  cmplx &operator-=(const cmplx &c);
  cmplx &Conjugate() {
    _im *= -1;
    return *this;
  }
  double _re;
  double _im;
};

inline cmplx &cmplx::operator*=(const double &d) {
  _re *= d;
  _im *= d;
  return *this;
}

inline cmplx &cmplx::operator+=(const cmplx &c) {
  _re += c._re;
  _im += c._im;
  return *this;
}

inline cmplx &cmplx::operator-=(const cmplx &c) {
  _re -= c._re;
  _im -= c._im;
  return *this;
}

inline cmplx operator+(const cmplx &cl, const cmplx &cr) {
  return (cmplx(cl) += cr);
}

inline cmplx operator-(const cmplx &cl, const cmplx &cr) {
  return (cmplx(cl) -= cr);
}

inline cmplx operator*(const cmplx &cl, const double &d) {
  return (cmplx(cl) *= d);
}

inline cmplx operator*(const double &d, const cmplx &cr) {
  return (cmplx(cr) *= d);
}

inline cmplx operator*(const cmplx &cl, const cmplx &cr) {
  return cmplx(cl.Re() * cr.Re() - cl.Im() * cr.Im(),
               cl.Re() * cr.Im() + cl.Im() * cr.Re());
}

template <typename NrTyp = double>
struct triple {
  triple() { ; }
  triple(NrTyp pp, NrTyp pu, NrTyp uu) : _pp(pp), _pu(pu), _uu(uu) { ; }
  triple(const triple<NrTyp> &t) : _pp(t._pp), _pu(t._pu), _uu(t._uu) { ; }
  triple &operator=(const triple &t) = default; 
  triple &operator*=(const double &d);
  triple &operator+=(const triple &t);
  triple &operator-=(const triple &t);
  NrTyp Sum() { return _pp + _pu + _uu; }
  NrTyp _pp;
  NrTyp _pu;
  NrTyp _uu;
};

template <typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator*=(const double &d) {
  _pp *= d;
  _pu *= d;
  _uu *= d;
  return *this;
}

template <typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator+=(const triple<NrTyp> &t) {
  _pp += t._pp;
  _pu += t._pu;
  _uu += t._uu;
  return *this;
}

template <typename NrTyp>
inline triple<NrTyp> &triple<NrTyp>::operator-=(const triple<NrTyp> &t) {
  _pp -= t._pp;
  _pu -= t._pu;
  _uu -= t._uu;
  return *this;
}

template <typename NrTyp>
inline triple<NrTyp> operator+(const triple<NrTyp> &tl,
                               const triple<NrTyp> &tr) {
  return (triple<NrTyp>(tl) += tr);
}

template <typename NrTyp>
inline triple<NrTyp> operator-(const triple<NrTyp> &tl,
                               const triple<NrTyp> &tr) {
  return (triple<NrTyp>(tl) -= tr);
}

template <typename NrTyp>
inline triple<NrTyp> operator*(const triple<NrTyp> &tl, const double &d) {
  return (triple<NrTyp>(tl) *= d);
}

template <typename NrTyp>
inline triple<NrTyp> operator*(const double &d, const triple<NrTyp> &tr) {
  return (triple<NrTyp>(tr) *= d);
}

template <typename NrTyp>
inline std::ostream &operator<<(std::ostream &out, const triple<NrTyp> &tr) {
  out << (boost::format("( %1$+1.7e %2$+1.7e %3$+1.7e )") % tr._pp % tr._pu %
          tr._uu);
  return out;
}

// To sort K-vectors via std::sort using a norm functor
template <class Norm, class V>
struct VectorSort {
  VectorSort() : _p(1e-40) { ; }
  VectorSort(double precision) : _p(precision) { ; }
  inline bool operator()(const V &v1, const V &v2);
  inline bool operator()(const V *v1, const V *v2);
  inline bool MatchDouble(double a, double b) {
    return ((a - b) * (a - b) < _p) ? true : false;
  }
  double _p;
  Norm _norm;
};

// Tschebyschow norm functor
struct MaxNorm {
  inline double operator()(const vec &v) { return v.cwiseAbs().maxCoeff(); }
};
// Euclidean norm functor
struct EucNorm {
  inline double operator()(const vec &v) { return v.norm(); }
};

// K-vector class (for grading purposes)
class KVector {
 public:
  KVector(vec k, double grade) : _k(k), _grade(grade), _has_sfactor(false) { ; }
  const vec &getK() const { return _k; }
  const double &getGrade() const { return _grade; }
  const double &getX() const { return _k(0); }
  const double &getY() const { return _k(1); }
  const double &getZ() const { return _k(2); }

  void setStructureFactor(EWD::cmplx &sfac) {
    _sfactor = sfac;
    _has_sfactor = true;
  }
  const EWD::cmplx &getStructureFactor() { return _sfactor; }
  double operator()(int a) { return _k(a); }

  double operator()(int a) const { return _k(a); }

 private:
  vec _k;
  double _grade;
  EWD::cmplx _sfactor;
  bool _has_sfactor;
};

// Specialized K-vector norm
struct KNorm {
  inline double operator()(const KVector &v) { return -v.getGrade(); }
  inline double operator()(const KVector *v) { return -v->getGrade(); }
};

template <class Norm, class V>
inline bool VectorSort<Norm, V>::operator()(const V &v1, const V &v2) {
  bool smaller = false;
  // LEVEL 1: MAGNITUDE
  double V1 = _norm(v1);
  double V2 = _norm(v2);
  if (MatchDouble(V1, V2)) {
    // LEVEL 2: X
    double X1 = v1(0);
    double X2 = v2(0);
    if (MatchDouble(X1, X2)) {
      // LEVEL 3: Y
      double Y1 = v1(1);
      double Y2 = v2(1);
      if (MatchDouble(Y1, Y2)) {
        // LEVEL 4: Z
        double Z1 = v1(2);
        double Z2 = v2(2);
        if (MatchDouble(Z1, Z2))
          smaller = true;
        else
          smaller = (Z1 < Z2) ? true : false;
      } else
        smaller = (Y1 < Y2) ? true : false;
    } else
      smaller = (X1 < X2) ? true : false;
  } else
    smaller = (V1 < V2) ? true : false;
  return smaller;
}

template <class Norm, class V>
inline bool VectorSort<Norm, V>::operator()(const V *v1, const V *v2) {
  bool smaller = false;
  // LEVEL 1: MAGNITUDE
  double V1 = _norm(v1);
  double V2 = _norm(v2);
  if (MatchDouble(V1, V2)) {
    // LEVEL 2: X
    double X1 = v1->getX();
    double X2 = v2->getX();
    if (MatchDouble(X1, X2)) {
      // LEVEL 3: Y
      double Y1 = v1->getY();
      double Y2 = v2->getY();
      if (MatchDouble(Y1, Y2)) {
        // LEVEL 4: Z
        double Z1 = v1->getZ();
        double Z2 = v2->getZ();
        if (MatchDouble(Z1, Z2))
          smaller = true;
        else
          smaller = (Z1 < Z2) ? true : false;
      } else
        smaller = (Y1 < Y2) ? true : false;
    } else
      smaller = (X1 < X2) ? true : false;
  } else
    smaller = (V1 < V2) ? true : false;
  return smaller;
}

}  // namespace EWD
}  // namespace xtp
}  // namespace votca

#endif
