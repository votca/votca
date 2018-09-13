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
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_DMASPACE_H
#define VOTCA_XTP_DMASPACE_H

#include <votca/tools/vec.h>
#include <votca/xtp/ewdspace.h>

#include <cmath>
#include <iostream>
#include <exception>

#include <boost/math/special_functions/binomial.hpp>

namespace votca {
namespace xtp {
namespace DMA {

// TODO Maybe move cmplx to separate header - or use C++ <complex> class
typedef votca::xtp::EWD::cmplx cmplx;
    
// BINOMIAL COEFFICIENT
inline double binom(int upper, int lower) {
    return boost::math::binomial_coefficient<double>(upper, lower);
}

// REGULAR SPHERICAL HARMONICS
// o NOTE Currently only up to rank l=2 (<> quadrupole)
// o TODO For higher ranks, use recursive construction
// o Notation key: C1_1 => l=1 m=-1 and C1x1 => l=1 m=+1
class RegularSphericalHarmonics
{
public:
    RegularSphericalHarmonics(votca::tools::vec r) 
       : _x(r.getX()), _y(r.getY()), _z(r.getZ()), _r2(r*r) {}    
   ~RegularSphericalHarmonics() {}
    
    // Real spherical harmonics
    inline double R00()  { return 1; }    
    inline double R10()  { return _z; }
    inline double R11c() { return _x; }
    inline double R11s() { return _y; }    
    inline double R20()  { return 0.5*(3*_z*_z - _r2); }
    inline double R21c() { return std::sqrt(3)*_x*_z; }
    inline double R21s() { return std::sqrt(3)*_y*_z; }
    inline double R22c() { return 0.5*std::sqrt(3)*(_x*_x - _y*_y); }
    inline double R22s() { return std::sqrt(3)*_x*_y; }
   
    // Complex spherical harmonics
    inline cmplx C00()   { return cmplx(R00(), 0.0); }    
    inline cmplx C1_1()  { return std::sqrt(0.5)*cmplx(+1*R11c(), -1*R11s()); }
    inline cmplx C10()   { return cmplx(R10(), 0.0); }
    inline cmplx C1x1()  { return std::sqrt(0.5)*cmplx(-1*R11c(), -1*R11s()); }    
    inline cmplx C2_2()  { return std::sqrt(0.5)*cmplx(+1*R22c(), -1*R22s()); }
    inline cmplx C2_1()  { return std::sqrt(0.5)*cmplx(+1*R21c(), -1*R21s()); }
    inline cmplx C20()   { return cmplx(R20(), 0.0); }
    inline cmplx C2x1()  { return std::sqrt(0.5)*cmplx(-1*R21c(), -1*R21s()); }
    inline cmplx C2x2()  { return std::sqrt(0.5)*cmplx(+1*R22c(), +1*R22s()); }
    
private:
    double _x, _y, _z;
    double _r2;    
};


// SPHERICAL MOMENT CONVERTER - REAL TO COMPLEX
// o Real spherical moments ordered like this:    
//   lm = 00    10 11c 11s    20 21c 21s 22c 22s ...
// o Complex spherical moments ordered like this: 
//   lm = 0+0   1-1 1+0 1+1   2-2 2-1 20 2+1 2+2 ...
class ComplexSphericalMoments
{
public:
    ComplexSphericalMoments(const std::vector<double> &Qlm) : _Qlm(Qlm) 
    { 
      if(_Qlm.size()!=9){
        throw std::invalid_argument("ComplexSphericalMoments must take a 9 value vector");
      }
    }
   ~ComplexSphericalMoments() {}
   
    // Real spherical moments (from constructor input)
    const double &Q00()  { return _Qlm[0]; }    
    const double &Q10()  { return _Qlm[1]; }
    const double &Q11c() { return _Qlm[2]; }
    const double &Q11s() { return _Qlm[3]; }    
    const double &Q20()  { return _Qlm[4]; }
    const double &Q21c() { return _Qlm[5]; }
    const double &Q21s() { return _Qlm[6]; }
    const double &Q22c() { return _Qlm[7]; }
    const double &Q22s() { return _Qlm[8]; }
   
    // Complex spherical moments
    inline cmplx X00()   { return cmplx(Q00(), 0.0); }    
    inline cmplx X1_1()  { return std::sqrt(0.5)*cmplx(+1*Q11c(), -1*Q11s()); }
    inline cmplx X10()   { return cmplx(Q10(), 0.0); }
    inline cmplx X1x1()  { return std::sqrt(0.5)*cmplx(-1*Q11c(), -1*Q11s()); }    
    inline cmplx X2_2()  { return std::sqrt(0.5)*cmplx(+1*Q22c(), -1*Q22s()); }
    inline cmplx X2_1()  { return std::sqrt(0.5)*cmplx(+1*Q21c(), -1*Q21s()); }
    inline cmplx X20()   { return cmplx(Q20(), 0.0); }
    inline cmplx X2x1()  { return std::sqrt(0.5)*cmplx(-1*Q21c(), -1*Q21s()); }
    inline cmplx X2x2()  { return std::sqrt(0.5)*cmplx(+1*Q22c(), +1*Q22s()); }
    
    void PrintReal() {
        std::cout << std::endl << std::scientific << "Q00 " << Q00();
        std::cout << std::endl << std::scientific << "Q10 " << Q10();
        std::cout << std::endl << std::scientific << "    " << Q11c();
        std::cout << std::endl << std::scientific << "    " << Q11s();
        std::cout << std::endl << std::scientific << "Q20 " << Q20();
        std::cout << std::endl << std::scientific << "    " << Q21c();
        std::cout << std::endl << std::scientific << "    " << Q21s();
        std::cout << std::endl << std::scientific << "    " << Q22c();
        std::cout << std::endl << std::scientific << "    " << Q22s();
    }
    
    std::vector<cmplx> ToVector() {
        std::vector<cmplx> Xlm;
        Xlm.push_back(X00());
        Xlm.push_back(X1_1());
        Xlm.push_back(X10());
        Xlm.push_back(X1x1());
        Xlm.push_back(X2_2());
        Xlm.push_back(X2_1());
        Xlm.push_back(X20());
        Xlm.push_back(X2x1());
        Xlm.push_back(X2x2());
        return Xlm;
    }
    
private:
    const std::vector<double> &_Qlm;
};


// SPHERICAL MOMENT CONVERTER - COMPLEX TO REAL
// o Real spherical moments ordered like this:    
//   lm = 00    10 11c 11s    20 21c 21s 22c 22s ...
// o Complex spherical moments ordered like this: 
//   lm = 0+0   1-1 1+0 1+1   2-2 2-1 20 2+1 2+2 ...
class RealSphericalMoments
{
public:
    RealSphericalMoments(const std::vector<cmplx> &Xlm) : _Xlm(Xlm) 
    { 
      if(_Xlm.size()!=9){
        throw std::invalid_argument("RealsphericalMoments class must take a 9 value vector");
      }
    }
   ~RealSphericalMoments() {}
   
    // Real spherical moments
    inline double Q00()  { return X00().Re(); }    
    inline double Q10()  { return X10().Re(); }
    inline double Q11c() { 
      cmplx q = std::sqrt(0.5) * cmplx(1,0) * (X1_1() - X1x1()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }
    inline double Q11s() { 
      cmplx q = sqrt(0.5) * cmplx(0,1) * (X1_1() + X1x1()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }    
    inline double Q20()  { return X20().Re(); }
    inline double Q21c() { 
      cmplx q = sqrt(0.5) * cmplx(1,0) * (X2_1() - X2x1()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }
    inline double Q21s() { 
      cmplx q = sqrt(0.5) * cmplx(0,1) * (X2_1() + X2x1()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }
    inline double Q22c() { 
      cmplx q = sqrt(0.5) * cmplx(1,0) * (X2_2() + X2x2()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }
    inline double Q22s() { 
      cmplx q = sqrt(0.5) * cmplx(0,1) * (X2_2() - X2x2()); 
      assert(q.Im() < 1e-8); 
      return q.Re(); 
    }

    // Complex spherical moments (from constructor input)
    const cmplx &X00()   { return _Xlm[0]; }    
    const cmplx &X1_1()  { return _Xlm[1]; }
    const cmplx &X10()   { return _Xlm[2]; }
    const cmplx &X1x1()  { return _Xlm[3]; }    
    const cmplx &X2_2()  { return _Xlm[4]; }
    const cmplx &X2_1()  { return _Xlm[5]; }
    const cmplx &X20()   { return _Xlm[6]; }
    const cmplx &X2x1()  { return _Xlm[7]; }
    const cmplx &X2x2()  { return _Xlm[8]; }
    
    void AddToVector(std::vector<double> &base) {        
      if(base.size()!=9){
        throw std::invalid_argument("AddToVector function must take a 9 value vector");
      }
      base[0] += Q00();
      base[1] += Q10();
      base[2] += Q11c();
      base[3] += Q11s();
      base[4] += Q20();
      base[5] += Q21c();
      base[6] += Q21s();
      base[7] += Q22c();
      base[8] += Q22s();
    }
    
    std::vector<double> ToVector() {
        std::vector<double> Qlm;
        Qlm.push_back(Q00());
        Qlm.push_back(Q10());
        Qlm.push_back(Q11c());
        Qlm.push_back(Q11s());
        Qlm.push_back(Q20());
        Qlm.push_back(Q21c());
        Qlm.push_back(Q21s());
        Qlm.push_back(Q22c());
        Qlm.push_back(Q22s());
        return Qlm;
    }
    
    void PrintReal() {
        std::cout << std::endl << std::scientific << "Q00 " << Q00();
        std::cout << std::endl << std::scientific << "Q10 " << Q10();
        std::cout << std::endl << std::scientific << "    " << Q11c();
        std::cout << std::endl << std::scientific << "    " << Q11s();
        std::cout << std::endl << std::scientific << "Q20 " << Q20();
        std::cout << std::endl << std::scientific << "    " << Q21c();
        std::cout << std::endl << std::scientific << "    " << Q21s();
        std::cout << std::endl << std::scientific << "    " << Q22c();
        std::cout << std::endl << std::scientific << "    " << Q22s();
    }
    
private:
    const std::vector<cmplx> &_Xlm;
};


class MomentShift
{
public:
    MomentShift() {}
   ~MomentShift() {}
    
    std::vector<cmplx> Shift(ComplexSphericalMoments &Xlm, 
        RegularSphericalHarmonics &Clm) {
        
        std::vector<cmplx> Xlm_shifted;
        
        // CHARGE
        cmplx X00 = 
                SqrtBinBin(0,0,0,0)   * Xlm.X00()*Clm.C00();  // L=0 M=0  l=0  m=0
        Xlm_shifted.push_back(X00);
        
        // DIPOLE
        cmplx X1_1 = 
                SqrtBinBin(1,-1,0,0)  * Xlm.X00()*Clm.C1_1()  // L=1 M=-1 l=0  m=0
              + SqrtBinBin(1,-1,1,-1) * Xlm.X1_1()*Clm.C00(); // L=1 M=-1 l=1  m=-1
        cmplx X10 =
                SqrtBinBin(1,0,0,0)   * Xlm.X00()*Clm.C10()   // L=1 M=0  l=0  m=0
              + SqrtBinBin(1,0,1,0)   * Xlm.X10()*Clm.C00();  // L=1 M=0  l=1  m=0        
        cmplx X1x1 =
                SqrtBinBin(1,1,0,0)   * Xlm.X00()*Clm.C1x1()  // L=1 M=+1 l=0  m=0
              + SqrtBinBin(1,1,1,1)   * Xlm.X1x1()*Clm.C00(); // L=1 M=+1 l=+1 m=+1
        Xlm_shifted.push_back(X1_1);
        Xlm_shifted.push_back(X10);
        Xlm_shifted.push_back(X1x1);
        
        // QUADRUPOLE
        cmplx X2_2 =
                SqrtBinBin(2,-2,0,0)  * Xlm.X00()*Clm.C2_2()  // L=2 M=-2 l=0  m=0
              + SqrtBinBin(2,-2,1,-1) * Xlm.X1_1()*Clm.C1_1() // L=2 M=-2 l=1  m=-1
              + SqrtBinBin(2,-2,2,-2) * Xlm.X2_2()*Clm.C00(); // L=2 M=-2 l=2  m=-2
        cmplx X2_1 =
                SqrtBinBin(2,-1,0,0)  * Xlm.X00()*Clm.C2_1()  // ...
              + SqrtBinBin(2,-1,1,0)  * Xlm.X10()*Clm.C1_1()
              + SqrtBinBin(2,-1,1,-1) * Xlm.X1_1()*Clm.C10()
              + SqrtBinBin(2,-1,2,-1) * Xlm.X2_1()*Clm.C00();
        cmplx X20 =
                SqrtBinBin(2,0,0,0)   * Xlm.X00()*Clm.C20()
              + SqrtBinBin(2,0,1,0)   * Xlm.X10()*Clm.C10()
              + SqrtBinBin(2,0,1,-1)  * Xlm.X1_1()*Clm.C1x1()
              + SqrtBinBin(2,0,1,+1)  * Xlm.X1x1()*Clm.C1_1()
              + SqrtBinBin(2,0,2,0)   * Xlm.X20()*Clm.C00();
        cmplx X2x1 =
                SqrtBinBin(2,1,0,0)   * Xlm.X00()*Clm.C2x1()
              + SqrtBinBin(2,1,1,0)   * Xlm.X10()*Clm.C1x1()
              + SqrtBinBin(2,1,1,1)   * Xlm.X1x1()*Clm.C10()
              + SqrtBinBin(2,1,2,1)   * Xlm.X2x1()*Clm.C00();
        cmplx X2x2 =
                SqrtBinBin(2,2,0,0)   * Xlm.X00()*Clm.C2x2()
              + SqrtBinBin(2,2,1,1)   * Xlm.X1x1()*Clm.C1x1()
              + SqrtBinBin(2,2,2,2)   * Xlm.X2x2()*Clm.C00();
        Xlm_shifted.push_back(X2_2);
        Xlm_shifted.push_back(X2_1);
        Xlm_shifted.push_back(X20);
        Xlm_shifted.push_back(X2x1);
        Xlm_shifted.push_back(X2x2);
        
        return Xlm_shifted;
    }
    
    double SqrtBinBin(int L, int M, int l, int m) {
        return std::sqrt( binom(L+M, l+m) * binom(L-M, l-m) );
    }
    
};

}}}

#endif // VOTCA_XTP_DMASPACE_H
