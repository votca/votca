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

#ifndef VOTCA_XTP_EWDINTERACTOR_H
#define	VOTCA_XTP_EWDINTERACTOR_H

#include <cmath>
#include <votca/tools/vec.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/polarseg.h>
#include <votca/xtp/ewdspace.h>

namespace votca { namespace xtp {

/// UNITS IN INPUT FILES
/// ... ... Positions as required by format
/// ... ... Multipole moment of rank k in e(a0)**k
/// ... ... Dipole polarizability in A³ (Angstrom cubed)

/// UNITS USED INTERNALLY
/// ... ... Positions in nm
/// ... ... Multipole moment of rank k in e(nm)**k
/// ... ... Dipole polarizability in nm³

/// CONVERSION FACTORS
/// ... ... Electric field (N/C) = Field (int)  * 1/4PiEps0(SI) * e * 1e+18
/// ... ... Energy (eV)          = Energy (int) * 1/4PiEps0(SI) * e * 1e+09
/// ... ... Potential(V)         = Pot. (int)   * 1/4PiEps0(SI) * e * 1e+09    
    
class EwdInteractor
{
public:

    EwdInteractor(double alpha, double tsharp) {
        a1 = alpha;
        a2 = a1*a1;
        a3 = a1*a2;
        a4 = a1*a3;
        a5 = a1*a4;
        a6 = a1*a5;
        a7 = a1*a6;
        
        ta1 = tsharp;
        ta2 = ta1*ta1;
        ta3 = ta1*ta2;
    };
    
    EwdInteractor() {};
   ~EwdInteractor() {};
   
    //static const double int2eV  = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    //static const double int2V_m = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-18;
    //static const double rSqrtPi = 0.564189583547756279280349644977832;    
    
    // Thole damping functions
    inline double L3() { return 1 - exp( -ta1*tu3); }
    inline double L5() { return 1 - (1 + ta1*tu3) * exp( -ta1*tu3); }
    inline double L7() { return 1 - (1 + ta1*tu3 + 0.6*ta2*tu3*tu3) * exp( -ta1*tu3); }
    inline double L9() { return 1 - (1 + ta1*tu3 + (18*ta2*tu3*tu3 + 9*ta3*tu3*tu3*tu3)/35) * exp( -ta1*tu3); }
    
    
    // ============================= REAL SPACE ============================= //    
    inline void ApplyBias(APolarSite &p1, APolarSite &p2);
    inline void ApplyBiasPolar(APolarSite &p1, APolarSite &p2);
    inline void ApplyBias(APolarSite &p1, APolarSite &p2, vec &s);
    inline void ApplyBiasPolar(APolarSite &p1, APolarSite &p2, vec &s);
    
    // =============================== //
    // ISOTROPIC & ANISOTROPIC KERNELS //
    // =============================== //
    
    // Make sure to set R1, R2, ... and rR1, rR2, ... before using {gB0, ...}
    inline void UpdateAllBls();
    inline double gB0() { return erfc(a1*R1)*rR1; }
    inline double gB1() { return rR2*(   gB0()  +  2*a1*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gB2() { return rR2*( 3*gB1()  +  4*a3*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gB3() { return rR2*( 5*gB2()  +  8*a5*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gB4() { return rR2*( 7*gB3()  + 16*a7*EWD::rSqrtPi * exp(-a2*R2)); }    
    // Make sure to set R1, R2, ... and rR1, rR2, ... before using {gC0, ...}
    inline void UpdateAllCls();
    inline double gC0() { return erf(a1*R1)*rR1; }
    inline double gC1() { return rR2*(   gC0()  -  2*a1*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gC2() { return rR2*( 3*gC1()  -  4*a3*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gC3() { return rR2*( 5*gC2()  -  8*a5*EWD::rSqrtPi * exp(-a2*R2)); }
    inline double gC4() { return rR2*( 7*gC3()  - 16*a7*EWD::rSqrtPi * exp(-a2*R2)); }    
    // Make sure to set R1, R2, ... and rxx, rxy, ... before using {gG0, ...}
    inline void UpdateAllGls(APolarSite &p1, APolarSite &p2);
    inline double gG0(APolarSite &p1, APolarSite &p2);
    inline double gG1(APolarSite &p1, APolarSite &p2);
    inline double gG2(APolarSite &p1, APolarSite &p2);
    inline double gG3(APolarSite &p1, APolarSite &p2);
    inline double gG4(APolarSite &p1, APolarSite &p2);
    
    // ====== //
    // FIELDS //
    // ====== //    
    
    // Real-space term contribution P1 <> P2
    inline double FPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    inline double FP12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    inline double FU12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    inline double FPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s);
    inline double FP12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s);
    inline double FU12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s);
    // Reciprocal-space double-counting correction term P1 <> P2
    inline double FPU12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    inline double FP12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    inline double FU12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    // Reciprocal-space K=0 shape correction term P1 <> P2
    void FPU12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);
    void FP12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);
    void FU12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);    
    
    // ========== //
    // POTENTIALS //
    // ========== //
    
    // Real-space term contribution P1 <> P2
    inline double PhiPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    inline double PhiP12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    inline double PhiU12_ERFC_At_By(APolarSite &p1, APolarSite &p2);
    // Reciprocal-space self-interaction term P1 <> P2
    inline double PhiPU12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    inline double PhiP12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    inline double PhiU12_ERF_At_By(APolarSite &p1, APolarSite &p2);
    // Reciprocal-space K=0 shape correction term P1 <> P2
    void PhiPU12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);
    void PhiP12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);
    void PhiU12_ShapeField_At_By(std::vector<PolarSeg*> &at, std::vector<PolarSeg*> &by, std::string shape, double volume);
    
    // ======== //
    // ENERGIES //
    // ======== //
    
    // Real-space term contribution P1 <> P2
    inline EWD::triple<double> U12_ERFC(APolarSite &p1, APolarSite &p2);    
    // Reciprocal-space double-counting correction term P1 <> P2
    inline EWD::triple<double> U12_ERF(APolarSite &p1, APolarSite &p2);    
    // Reciprocal-space K=0 shape correction term P1 <> P2
    EWD::triple<double> U12_ShapeTerm(std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, std::string shape, double volume, Logger *log=NULL);
    
    // ========================== RECIPROCAL SPACE ========================== //    
    inline void   ApplyBiasK(const vec &k);
    inline double Ark2Expk2(const vec &k);
    
    inline void PUApplyBiasK(APolarSite &p);
    inline void PApplyBiasK(APolarSite &p);
    inline void UApplyBiasK(APolarSite &p);
    
    // Structure amplitudes //
    EWD::cmplx PUStructureAmplitude(std::vector<PolarSeg*> &s);
    EWD::cmplx PStructureAmplitude(std::vector<PolarSeg*> &s);
    EWD::cmplx UStructureAmplitude(std::vector<PolarSeg*> &s);    
    EWD::cmplx PUStructureAmplitude(std::vector<PolarSeg*> &s, const vec &k);
    EWD::cmplx PStructureAmplitude(std::vector<PolarSeg*> &s, const vec &k);
    EWD::cmplx UStructureAmplitude(std::vector<PolarSeg*> &s, const vec &k);
    
    // Fields
    EWD::cmplx FPU12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);
    EWD::cmplx FP12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);
    EWD::cmplx FU12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);    

    EWD::cmplx FP12_At_ByS2(const vec &k, std::vector<PolarSeg*> &s1, const EWD::cmplx &S2, double &rV);
    EWD::cmplx FU12_At_ByS2(const vec &k, std::vector<PolarSeg*> &s1, const EWD::cmplx &S2, double &rV);
    
    // Potentials
    EWD::cmplx PhiPU12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);
    EWD::cmplx PhiP12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);
    EWD::cmplx PhiU12_AS1S2_At_By(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2, double &rV);
    
    // Energies
    EWD::triple<EWD::cmplx> AS1S2(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2);
    EWD::triple<EWD::cmplx> S1S2(const vec &k, std::vector<PolarSeg*> &s1, std::vector<PolarSeg*> &s2);
    
private:
    
    // Thole sharpness parameter & reduced interaction distance
    double ta1, ta2, ta3;
    double tu3 = 0;
    double l3 = 0, l5 = 0, l7 = 0, l9 = 0;     // (damps fields originating in indu. m'poles)
    double lp3 = 0, lp5 = 0, lp7 = 0, lp9 = 0; // (damps fields originating in perm. m'poles)
    
    // Ewald sharpness parameter powers
    double a1, a2, a3, a4, a5, a6, a7;
    
    // ============================= REAL SPACE ============================= //
    
    // Connection vector (1) <- (2), i.e. r12 = r1 - r2 == rab = ra - rb;
    vec r12 = 0.0;
    
    // Vector components rx = r12x, ...
    double rx = 0.0, ry = 0.0, rz = 0.0;
    
    // Matrix product rxx = r12x*r12x, rxy = r12x*r12y, ...
    double rxx = 0.0, rxy = 0.0, rxz = 0.0, ryy = 0.0, ryz = 0.0, rzz = 0.0;
    
    // Real-space distance powers
    double R1 = 0.0, R2 = 0.0, R3 = 0.0;//, R4, R5;
    
    // Real-space inverse distance powers
    double rR1 = 0.0, rR2 = 0.0;//, rR3, rR4, rR5;
    
    // {G} function values
    double ppG0 = 0.0, ppG1 = 0.0, ppG2 = 0.0, ppG3 = 0.0, ppG4 = 0.0;
    double puG1 = 0.0, puG2 = 0.0, puG3 = 0.0;
    double uuG1 = 0.0, uuG2 = 0.0;
    
    // {Bl}, {Cl} function values
    double rSqrtPiExp = 0.0;
    double B0 = 0.0, B1 = 0.0, B2 = 0.0, B3 = 0.0, B4 = 0.0;
    double C0 = 0.0, C1 = 0.0, C2 = 0.0, C3 = 0.0, C4 = 0.0;
    
    
    // ========================== RECIPROCAL SPACE ========================== //
    
    // k-space vector
    vec k12 = 0.0;
    double K = 0.0;
    double AK = 0.0;
    
    // Vector components kx, ...
    double kx = 0.0, ky = 0.0, kz = 0.0;
    
    // Matrix product kxx = kx*kx, kxy = kx*ky, ...
    double kxx = 0.0, kxy = 0.0, kxz = 0.0, kyy = 0.0, kyz = 0.0, kzz = 0.0;
    
    // cos(k*r), sin(k*r), µ * k, Q : K
    double dk = 0.0;             // <- µ(perm.) * k
    double u_dk = 0.0;           // <- µ(indu.) * k
    double Qk = 0.0;
    double kr = 0.0;
    double coskr = 0.0;
    double sinkr = 0.0;    
    double re_s = 0.0, im_s = 0.0;     // <- perm. moments
    double u_re_s = 0.0, u_im_s = 0.0; // <- indu. moments
};

// =============================== REAL SPACE =============================== //

inline void EwdInteractor::ApplyBias(APolarSite& p1, APolarSite& p2) {
    
    r12 = p1.getPos() - p2.getPos();
    
    rx = r12.getX();
    ry = r12.getY();
    rz = r12.getZ();
    
    rxx = rx*rx;     rxy = rx*ry;     rxz = rx*rz;
    ryy = ry*ry;     ryz = ry*rz;
    rzz = rz*rz;
                                      
    R1 = votca::tools::abs(r12);
    R2 = R1*R1;
    //R3 = R1*R2;
    //R4 = R1*R3;
    //R5 = R1*R4;
    
    rR1 = 1./R1;
    rR2 = 1./R2;
    //rR3 = 1./R3;
    //rR4 = 1./R4;
    //rR5 = 1./R5;
    
    return;
}


inline void EwdInteractor::ApplyBiasPolar(APolarSite& p1, APolarSite& p2) {
    
    r12 = p1.getPos() - p2.getPos();
    
    rx = r12.getX();
    ry = r12.getY();
    rz = r12.getZ();
    
    rxx = rx*rx;     rxy = rx*ry;     rxz = rx*rz;
    ryy = ry*ry;     ryz = ry*rz;
    rzz = rz*rz;
                                      
    R1 = votca::tools::abs(r12);
    R2 = R1*R1;
    R3 = R1*R2;
    //R4 = R1*R3;
    //R5 = R1*R4;
    
    rR1 = 1./R1;
    rR2 = 1./R2;
    //rR3 = 1./R3;
    //rR4 = 1./R4;
    //rR5 = 1./R5;
    
    // Thole damping initialization
    tu3   = R3 / sqrt(
        1./3.*(p1.Pxx*p2.Pxx + p1.Pxy*p2.Pxy + p1.Pxz*p2.Pxz
             + p1.Pxy*p2.Pxy + p1.Pyy*p2.Pyy + p1.Pyz*p2.Pyz
             + p1.Pxz*p2.Pxz + p1.Pyz*p2.Pyz + p1.Pzz*p2.Pzz) );
    
    // Thole damping functions
    if (ta1*tu3 < 40) {
        l3 = L3(); 
        l5 = L5();
        l7 = L7();
        l9 = L9();
        // cg <> cg interactions : do not damp PU terms (field & energy)
        if (false && p1._resolution == APolarSite::coarsegrained
         && p2._resolution == APolarSite::coarsegrained) {
                lp3 = lp5 = lp7 = lp9 = 1.;
        }
        else {
            lp3 = lp5 = lp7 = lp9 = 1.;
//            lp3 = l3;
//            lp5 = l5;
//            lp7 = l7;
//            lp9 = l9;
        }
    }
    else {
        l3 = l5 = l7 = l9 = 1.;
        lp3 = lp5 = lp7 = lp9 = 1.;
    }
    return;
}


inline void EwdInteractor::ApplyBias(APolarSite& p1, APolarSite& p2, vec &s) {
    // The shift vector s should be given as the following difference:
    //     s = dr(1->2,pbc) - dr(1->2),
    // where 'pbc' is short-hand for 'periodic-boundary corrected'. 
    // The effective position of particle 2 is hence r(2) + s.
    r12 = p1.getPos() - p2.getPos() - s;
    
    rx = r12.getX();
    ry = r12.getY();
    rz = r12.getZ();
    
    rxx = rx*rx;     rxy = rx*ry;     rxz = rx*rz;
    ryy = ry*ry;     ryz = ry*rz;
    rzz = rz*rz;
                                      
    R1 = votca::tools::abs(r12);
    R2 = R1*R1;
    //R3 = R1*R2;
    //R4 = R1*R3;
    //R5 = R1*R4;
    
    rR1 = 1./R1;
    rR2 = 1./R2;
    //rR3 = 1./R3;
    //rR4 = 1./R4;
    //rR5 = 1./R5;
    
    return;
}


inline void EwdInteractor::ApplyBiasPolar(APolarSite& p1, APolarSite& p2, vec &s) {
    // The shift vector s should be given as the following difference:
    //     s = dr(1->2,pbc) - dr(1->2),
    // where 'pbc' is short-hand for 'periodic-boundary corrected'. 
    // The effective position of particle 2 is hence r(2) + s.
    r12 = p1.getPos() - p2.getPos() - s;
    
    rx = r12.getX();
    ry = r12.getY();
    rz = r12.getZ();
    
    rxx = rx*rx;     rxy = rx*ry;     rxz = rx*rz;
    ryy = ry*ry;     ryz = ry*rz;
    rzz = rz*rz;
                                      
    R1 = votca::tools::abs(r12);
    R2 = R1*R1;
    R3 = R1*R2;
    //R4 = R1*R3;
    //R5 = R1*R4;
    
    rR1 = 1./R1;
    rR2 = 1./R2;
    //rR3 = 1./R3;
    //rR4 = 1./R4;
    //rR5 = 1./R5;
    
    // Thole damping initialization
    tu3   = R3 / sqrt(
        1./3.*(p1.Pxx*p2.Pxx + p1.Pxy*p2.Pxy + p1.Pxz*p2.Pxz
             + p1.Pxy*p2.Pxy + p1.Pyy*p2.Pyy + p1.Pyz*p2.Pyz
             + p1.Pxz*p2.Pxz + p1.Pyz*p2.Pyz + p1.Pzz*p2.Pzz) );
    
    // Thole damping functions
    if (ta1*tu3 < 40) {
        l3 = L3(); 
        l5 = L5();
        l7 = L7();
        l9 = L9();
        // cg <> cg interactions : do not damp PU terms (field & energy)
        if (false && p1._resolution == APolarSite::coarsegrained
         && p2._resolution == APolarSite::coarsegrained) {
                lp3 = lp5 = lp7 = lp9 = 1.;
        }
        else {
            lp3 = lp5 = lp7 = lp9 = 1.;
//            lp3 = l3;
//            lp5 = l5;
//            lp7 = l7;
//            lp9 = l9;
        }
    }
    else {
        l3 = l5 = l7 = l9 = 1.;
        lp3 = lp5 = lp7 = lp9 = 1.;
    }
    return;
}


inline void EwdInteractor::UpdateAllBls() {
    
    rSqrtPiExp = EWD::rSqrtPi * exp(-a2*R2);
    
    B0 = erfc(a1*R1)*rR1;    
    B1 = rR2*(   B0  +  2*a1*rSqrtPiExp);
    B2 = rR2*( 3*B1  +  4*a3*rSqrtPiExp);
    B3 = rR2*( 5*B2  +  8*a5*rSqrtPiExp);
    B4 = rR2*( 7*B3  + 16*a7*rSqrtPiExp);
    
//    B0 = gB0();
//    B1 = gB1();
//    B2 = gB2();
//    B3 = gB3();
//    B4 = gB4();
    
    return;
}


inline void EwdInteractor::UpdateAllCls() {
    
    double rSqrtPiExp = EWD::rSqrtPi * exp(-a2*R2);
    
    C0 = erf(a1*R1)*rR1;    
    C1 = rR2*(   C0  -  2*a1*rSqrtPiExp);
    C2 = rR2*( 3*C1  -  4*a3*rSqrtPiExp);
    C3 = rR2*( 5*C2  -  8*a5*rSqrtPiExp);
    C4 = rR2*( 7*C3  - 16*a7*rSqrtPiExp);

//    C0 = gC0();
//    C1 = gC1();
//    C2 = gC2();
//    C3 = gC3();
//    C4 = gC4();
    
    return;
}


inline void EwdInteractor::UpdateAllGls(APolarSite& p1, APolarSite& p2) {
    
    ppG0 = ppG1 = ppG2 = ppG3 = ppG4 = 0.0;
    puG1 = puG2 = puG3 = 0.0;
    uuG1 = uuG2 = 0.0;
    
    // Dot product µ * r
    double mu1_r = 0.0;
    double mu2_r = 0.0;
    
    // Dot product µ(ind) * r
    double u_mu1_r = 0.0;
    double u_mu2_r = 0.0;
    
    // Dot product Q * r
    double Q1_rx, Q1_ry, Q1_rz = 0.0;
    double Q2_rx, Q2_ry, Q2_rz = 0.0;
    
    // Dyadic product Q : R
    double Q1__R = 0.0;
    double Q2__R = 0.0;
    
    u_mu1_r = (p1.U1x*rx + p1.U1y*ry + p1.U1z*rz);
    if (p1._rank > 0) {        
        mu1_r = (p1.Q1x*rx + p1.Q1y*ry + p1.Q1z*rz);        
        if (p1._rank > 1) {
            Q1__R = (p1.Qxx*rxx + 2*p1.Qxy*rxy + 2*p1.Qxz*rxz
                                +   p1.Qyy*ryy + 2*p1.Qyz*ryz
                                               +   p1.Qzz*rzz);
            Q1_rx = p1.Qxx*rx + p1.Qxy*ry + p1.Qxz*rz;
            Q1_ry = p1.Qxy*rx + p1.Qyy*ry + p1.Qyz*rz;
            Q1_rz = p1.Qxz*rx + p1.Qyz*ry + p1.Qzz*rz;
        }
    }
    
    u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (p2._rank > 0) {
        mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        if (p2._rank > 1) {
            Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                               +   p2.Qzz*rzz);
            Q2_rx = p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz;
            Q2_ry = p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz;
            Q2_rz = p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz;
        }
    }    
    
    // RANK(1) >= RANK(2)
    // 1 - charge, 2 -charge
    ppG0 += p1.Q00 * p2.Q00;
    
    puG1 += -p2.Q00 * u_mu1_r;    
    uuG1 += p1.U1x*p2.U1x + p1.U1y*p2.U1y + p1.U1z*p2.U1z;
    uuG2 += - u_mu1_r * u_mu2_r;
    
    // 1 - dipole, 2 - charge
    if (p1._rank > 0) {
        ppG1 += - p2.Q00 * mu1_r;
        puG1 += p1.Q1x*p2.U1x + p1.Q1y*p2.U1y + p1.Q1z*p2.U1z;
        puG2 += - mu1_r * u_mu2_r;
        
        // 1 - dipole, 2 - dipole
        if (p2._rank > 0) {
            ppG1 += p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z;
            ppG2 += - mu1_r * mu2_r;
        }
        
        // 1 - quadrupole, 2 -charge
        if (p1._rank > 1) {
            ppG2 += p2.Q00 * Q1__R;
            
            // 1 - quadrupole, 2 - dipole
            if (p2._rank > 0) {
                ppG2 += - 2 * (p1.Qxx*p2.Q1x*rx + p1.Qxy*p2.Q1x*ry + p1.Qxz*p2.Q1x*rz
                             + p1.Qxy*p2.Q1y*rx + p1.Qyy*p2.Q1y*ry + p1.Qyz*p2.Q1y*rz
                             + p1.Qxz*p2.Q1z*rx + p1.Qyz*p2.Q1z*ry + p1.Qzz*p2.Q1z*rz);
                puG2 += - 2 * (p1.Qxx*p2.U1x*rx + p1.Qxy*p2.U1x*ry + p1.Qxz*p2.U1x*rz
                             + p1.Qxy*p2.U1y*rx + p1.Qyy*p2.U1y*ry + p1.Qyz*p2.U1y*rz
                             + p1.Qxz*p2.U1z*rx + p1.Qyz*p2.U1z*ry + p1.Qzz*p2.U1z*rz);
                ppG3 += mu2_r * Q1__R;
                puG3 += u_mu2_r * Q1__R;
                
                // 1 - quadrupole, 2 - quadrupole
                if (p2._rank > 1) {
                    ppG2 += 2 * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
                                               +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
                                                                 +   p1.Qzz*p2.Qzz);
                    ppG3 += -4 * (Q1_rx*Q2_rx + Q1_ry*Q2_ry + Q1_rz*Q2_rz);
                    ppG4 += Q1__R * Q2__R;
                }
            }            
        }        
    }
    
    
    // RANK(1) < RANK(2)    
    puG1 += + p1.Q00 * u_mu2_r;
    
    // 2 - dipole, 1 - charge
    if (p2._rank > 0) {
        ppG1 += + p1.Q00 * mu2_r;   
        puG1 += p1.U1x*p2.Q1x + p1.U1y*p2.Q1y + p1.U1z*p2.Q1z;
        puG2 += - u_mu1_r * mu2_r;
        
        // 2 - dipole, 1 - dipole
        // ... covered above.
        
        // 2 - quadrupole, 1 - charge
        if (p2._rank > 1) {
            ppG2 += p1.Q00 * Q2__R;
            
            // 2 - quadrupole, 1 - dipole
            if (p1._rank > 0) {
                ppG2 +=   2 * (p2.Qxx*p1.Q1x*rx + p2.Qxy*p1.Q1x*ry + p2.Qxz*p1.Q1x*rz
                             + p2.Qxy*p1.Q1y*rx + p2.Qyy*p1.Q1y*ry + p2.Qyz*p1.Q1y*rz
                             + p2.Qxz*p1.Q1z*rx + p2.Qyz*p1.Q1z*ry + p2.Qzz*p1.Q1z*rz);
                puG2 +=   2 * (p2.Qxx*p1.U1x*rx + p2.Qxy*p1.U1x*ry + p2.Qxz*p1.U1x*rz
                             + p2.Qxy*p1.U1y*rx + p2.Qyy*p1.U1y*ry + p2.Qyz*p1.U1y*rz
                             + p2.Qxz*p1.U1z*rx + p2.Qyz*p1.U1z*ry + p2.Qzz*p1.U1z*rz);
                ppG3 += - mu1_r * Q2__R;
                puG3 += - u_mu1_r * Q2__R;
                
                // 2 - quadrupole, 2 - quadrupole
                // ... covered above.
            }
        }
    }
    
    return;    
}


// =============================== REAL SPACE =============================== //
//                               ELECTRIC FIELD                               //

inline double EwdInteractor::FPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    // ATTENTION Increments p1->FP only, not p1->FU
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Charge
    fx += - p2.Q00*rx*lp3*B1;
    fy += - p2.Q00*ry*lp3*B1;
    fz += - p2.Q00*rz*lp3*B1;
    
    if (p2._rank > 0) {
        // Dipole
        fx += p2.Q1x*lp3*B1;
        fy += p2.Q1y*lp3*B1;
        fz += p2.Q1z*lp3*B1;
        
        double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        fx += - rx*mu2_r*lp5*B2;
        fy += - ry*mu2_r*lp5*B2;
        fz += - rz*mu2_r*lp5*B2;
        
        if (p2._rank > 1) {
            // Quadrupole
            fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * lp5*B2;
            fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * lp5*B2;
            fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * lp5*B2;
            
            double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                      +   p2.Qzz*rzz);            
            fx += - Q2__R*rx*lp7*B3;
            fy += - Q2__R*ry*lp7*B3;
            fz += - Q2__R*rz*lp7*B3;
        }
    }
    
    // Field generated by induced moments
    double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (ta1*tu3 < 40) {        
        fx += p2.U1x*l3*B1;
        fy += p2.U1y*l3*B1;
        fz += p2.U1z*l3*B1;        
        
        fx += - rx*u_mu2_r*l5*B2;
        fy += - ry*u_mu2_r*l5*B2;
        fz += - rz*u_mu2_r*l5*B2;
    }
    else {
        fx += p2.U1x*B1;
        fy += p2.U1y*B1;
        fz += p2.U1z*B1;        
        
        fx += - rx*u_mu2_r*B2;
        fy += - ry*u_mu2_r*B2;
        fz += - rz*u_mu2_r*B2;
    }
    
    // Increment fields
    p1.FPx += fx;
    p1.FPy += fy;
    p1.FPz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FP12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Charge
    fx += - p2.Q00*rx*lp3*B1;
    fy += - p2.Q00*ry*lp3*B1;
    fz += - p2.Q00*rz*lp3*B1;
    
    if (p2._rank > 0) {
        // Dipole
        fx += p2.Q1x*lp3*B1;
        fy += p2.Q1y*lp3*B1;
        fz += p2.Q1z*lp3*B1;
        
        double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        fx += - rx*mu2_r*lp5*B2;
        fy += - ry*mu2_r*lp5*B2;
        fz += - rz*mu2_r*lp5*B2;
        
        if (p2._rank > 1) {
            // Quadrupole
            fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * lp5*B2;
            fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * lp5*B2;
            fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * lp5*B2;
            
            double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                      +   p2.Qzz*rzz);            
            fx += - Q2__R*rx*lp7*B3;
            fy += - Q2__R*ry*lp7*B3;
            fz += - Q2__R*rz*lp7*B3;
        }
    }
    
    // Increment fields
    p1.FPx += fx;
    p1.FPy += fy;
    p1.FPz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FU12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Field generated by induced moments
    double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (ta1*tu3 < 40) {        
        fx += p2.U1x*l3*B1;
        fy += p2.U1y*l3*B1;
        fz += p2.U1z*l3*B1;        
        
        fx += - rx*u_mu2_r*l5*B2;
        fy += - ry*u_mu2_r*l5*B2;
        fz += - rz*u_mu2_r*l5*B2;
    }
    else {
        fx += p2.U1x*B1;
        fy += p2.U1y*B1;
        fz += p2.U1z*B1;        
        
        fx += - rx*u_mu2_r*B2;
        fy += - ry*u_mu2_r*B2;
        fz += - rz*u_mu2_r*B2;
    }
    
    // Increment fields
    p1.FUx += fx;
    p1.FUy += fy;
    p1.FUz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    // ATTENTION Increments p1->FP only, not p1->FU
    // The shift vector s should be given as the following difference:
    //     s = dr(1->2,pbc) - dr(1->2),
    // where 'pbc' is short-hand for 'periodic-boundary corrected'. 
    // The effective position of particle 2 is hence r(2) + s.
    ApplyBiasPolar(p1, p2, s);
    UpdateAllBls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Charge
    fx += - p2.Q00*rx*lp3*B1;
    fy += - p2.Q00*ry*lp3*B1;
    fz += - p2.Q00*rz*lp3*B1;
    
    if (p2._rank > 0) {
        // Dipole
        fx += p2.Q1x*lp3*B1;
        fy += p2.Q1y*lp3*B1;
        fz += p2.Q1z*lp3*B1;
        
        double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        fx += - rx*mu2_r*lp5*B2;
        fy += - ry*mu2_r*lp5*B2;
        fz += - rz*mu2_r*lp5*B2;
        
        if (p2._rank > 1) {
            // Quadrupole
            fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * lp5*B2;
            fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * lp5*B2;
            fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * lp5*B2;
            
            double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                      +   p2.Qzz*rzz);            
            fx += - Q2__R*rx*lp7*B3;
            fy += - Q2__R*ry*lp7*B3;
            fz += - Q2__R*rz*lp7*B3;
        }
    }
    
    // Field generated by induced moments
    double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (ta1*tu3 < 40) {
        fx += p2.U1x*l3*B1;
        fy += p2.U1y*l3*B1;
        fz += p2.U1z*l3*B1;        
        
        fx += - rx*u_mu2_r*l5*B2;
        fy += - ry*u_mu2_r*l5*B2;
        fz += - rz*u_mu2_r*l5*B2;
    }
    else {
        fx += p2.U1x*B1;
        fy += p2.U1y*B1;
        fz += p2.U1z*B1;        
        
        fx += - rx*u_mu2_r*B2;
        fy += - ry*u_mu2_r*B2;
        fz += - rz*u_mu2_r*B2;
    }
    
    // Increment fields
    p1.FPx += fx;
    p1.FPy += fy;
    p1.FPz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FP12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    // The shift vector s should be given as the following difference:
    //     s = dr(1->2,pbc) - dr(1->2),
    // where 'pbc' is short-hand for 'periodic-boundary corrected'. 
    // The effective position of particle 2 is hence r(2) + s.
    ApplyBiasPolar(p1, p2, s);
    UpdateAllBls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Charge
    fx += - p2.Q00*rx*lp3*B1;
    fy += - p2.Q00*ry*lp3*B1;
    fz += - p2.Q00*rz*lp3*B1;
    
    if (p2._rank > 0) {
        // Dipole
        fx += p2.Q1x*lp3*B1;
        fy += p2.Q1y*lp3*B1;
        fz += p2.Q1z*lp3*B1;
        
        double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        fx += - rx*mu2_r*lp5*B2;
        fy += - ry*mu2_r*lp5*B2;
        fz += - rz*mu2_r*lp5*B2;
        
        if (p2._rank > 1) {
            // Quadrupole
            fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * lp5*B2;
            fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * lp5*B2;
            fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * lp5*B2;
            
            double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                      +   p2.Qzz*rzz);            
            fx += - Q2__R*rx*lp7*B3;
            fy += - Q2__R*ry*lp7*B3;
            fz += - Q2__R*rz*lp7*B3;
        }
    }
    
    // Increment fields
    p1.FPx += fx;
    p1.FPy += fy;
    p1.FPz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FU12_ERFC_At_By(APolarSite &p1, APolarSite &p2, vec &s) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    // The shift vector s should be given as the following difference:
    //     s = dr(1->2,pbc) - dr(1->2),
    // where 'pbc' is short-hand for 'periodic-boundary corrected'. 
    // The effective position of particle 2 is hence r(2) + s.
    ApplyBiasPolar(p1, p2, s);
    UpdateAllBls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    // Field generated by induced moments
    double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (ta1*tu3 < 40) {
        fx += p2.U1x*l3*B1;
        fy += p2.U1y*l3*B1;
        fz += p2.U1z*l3*B1;        
        
        fx += - rx*u_mu2_r*l5*B2;
        fy += - ry*u_mu2_r*l5*B2;
        fz += - rz*u_mu2_r*l5*B2;
    }
    else {
        fx += p2.U1x*B1;
        fy += p2.U1y*B1;
        fz += p2.U1z*B1;        
        
        fx += - rx*u_mu2_r*B2;
        fy += - ry*u_mu2_r*B2;
        fz += - rz*u_mu2_r*B2;
    }
    
    // Increment fields
    p1.FUx += fx;
    p1.FUy += fy;
    p1.FUz += fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FPU12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    // ATTENTION Increments p1->FP only, not p1->FU
    ApplyBias(p1, p2);
    UpdateAllCls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    if (R1 < 1e-2) {
        if (p2._rank > 0) {
            fx += 4./3.*a3*EWD::rSqrtPi * p2.Q1x;
            fy += 4./3.*a3*EWD::rSqrtPi * p2.Q1y;
            fz += 4./3.*a3*EWD::rSqrtPi * p2.Q1z;
        }
        // Field generated by induced moments
        fx += 4./3.*a3*EWD::rSqrtPi * p2.U1x;
        fy += 4./3.*a3*EWD::rSqrtPi * p2.U1y;
        fz += 4./3.*a3*EWD::rSqrtPi * p2.U1z;
    }
    else {
        // Charge
        fx += - p2.Q00*rx*C1;
        fy += - p2.Q00*ry*C1;
        fz += - p2.Q00*rz*C1;    

        if (p2._rank > 0) {
            // Dipole
            fx += p2.Q1x*C1;
            fy += p2.Q1y*C1;
            fz += p2.Q1z*C1;

            double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
            fx += - rx*mu2_r*C2;
            fy += - ry*mu2_r*C2;
            fz += - rz*mu2_r*C2;

            if (p2._rank > 1) {
                // Quadrupole
                fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * C2;
                fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * C2;
                fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * C2;

                double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                           +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                          +   p2.Qzz*rzz);            
                fx += - Q2__R*rx*C3;
                fy += - Q2__R*ry*C3;
                fz += - Q2__R*rz*C3;
            }
        }
        // Field generated by induced moments
        double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
        
        fx += p2.U1x*C1;
        fy += p2.U1y*C1;
        fz += p2.U1z*C1;
        
        fx += - rx*u_mu2_r*C2;
        fy += - ry*u_mu2_r*C2;
        fz += - rz*u_mu2_r*C2;
    }
    // Increment fields. Note the (-): This is a compensation term.
    p1.FPx += - fx;
    p1.FPy += - fy;
    p1.FPz += - fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FP12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    ApplyBias(p1, p2);
    UpdateAllCls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    if (R1 < 1e-2) {
        if (p2._rank > 0) {
            fx += 4./3.*a3*EWD::rSqrtPi * p2.Q1x;
            fy += 4./3.*a3*EWD::rSqrtPi * p2.Q1y;
            fz += 4./3.*a3*EWD::rSqrtPi * p2.Q1z;
        }
    }
    else {
        // Charge
        fx += - p2.Q00*rx*C1;
        fy += - p2.Q00*ry*C1;
        fz += - p2.Q00*rz*C1;    

        if (p2._rank > 0) {
            // Dipole
            fx += p2.Q1x*C1;
            fy += p2.Q1y*C1;
            fz += p2.Q1z*C1;

            double mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
            fx += - rx*mu2_r*C2;
            fy += - ry*mu2_r*C2;
            fz += - rz*mu2_r*C2;

            if (p2._rank > 1) {
                // Quadrupole
                fx += 2 * (p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz) * C2;
                fy += 2 * (p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz) * C2;
                fz += 2 * (p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) * C2;

                double Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                           +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                          +   p2.Qzz*rzz);            
                fx += - Q2__R*rx*C3;
                fy += - Q2__R*ry*C3;
                fz += - Q2__R*rz*C3;
            }
        }
    }
    // Increment fields. Note the (-): This is a compensation term.
    p1.FPx += - fx;
    p1.FPy += - fy;
    p1.FPz += - fz;    
    return fx*fx + fy*fy + fz*fz;
}


inline double EwdInteractor::FU12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
    ApplyBias(p1, p2);
    UpdateAllCls();
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    if (R1 < 1e-2) {
        // Field generated by induced moments
        fx += 4./3.*a3*EWD::rSqrtPi * p2.U1x;
        fy += 4./3.*a3*EWD::rSqrtPi * p2.U1y;
        fz += 4./3.*a3*EWD::rSqrtPi * p2.U1z;
    }
    else {
        // Field generated by induced moments
        double u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
        
        fx += p2.U1x*C1;
        fy += p2.U1y*C1;
        fz += p2.U1z*C1;
        
        fx += - rx*u_mu2_r*C2;
        fy += - ry*u_mu2_r*C2;
        fz += - rz*u_mu2_r*C2;
    }
    // Increment fields. Note the (-): This is a compensation term.
    p1.FUx += - fx;
    p1.FUy += - fy;
    p1.FUz += - fz;    
    return fx*fx + fy*fy + fz*fz;
}


//inline double EwdInteractor::F12_XYSlab_At_By(APolarSite &p1, APolarSite &p2, 
//    double &TwoPi_V) {
//    assert(false); // This function is correct, but not necessarily efficient,
//                   // since an (external) double sum is not avoided
//    // NOTE Field points from (-) to (+) => XInductor, XInteractor compatible
//    ApplyBias(p1, p2);
//
//    double fz = 0.0;
//    
//    // Charge
//    fz += - 2*p2.Q00*rz;
//    if (p2._rank > 0) {
//        // Dipole
//        fz += 2*p2.Q1z;
//    }
//    
//    // Field generated by induced moments
//    fz += 2*p2.U1z;
//    
//    // Increment field, use prefactor 2*PI/V
//    p1.FPz += TwoPi_V*fz;    
//    return fz*fz;
//}




// =============================== REAL SPACE =============================== //
//                                 POTENTIALS                                 //


inline double EwdInteractor::PhiPU12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // ATTENTION Potentials are currently not damped
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    
    double phi_p = 0.0;
    double phi_u = 0.0;
    
    // Dot product µ * r
    double mu2_r = 0.0;
    // Dot product µ(ind) * r
    double u_mu2_r = 0.0;
    // Dyadic product Q : R
    double Q2__R = 0.0;
    
    u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    if (p2._rank > 0) {
        mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        if (p2._rank > 1) {
            Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                               +   p2.Qzz*rzz);
        }
    }
    
    phi_p = p2.Q00*B0 + mu2_r*B1 + Q2__R*B2;
    phi_u = u_mu2_r*B1;    
    p1.PhiP += phi_p;
    p1.PhiU += phi_u;    
    return phi_p + phi_u;
}


inline double EwdInteractor::PhiP12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // ATTENTION Potentials are currently not damped
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    
    double phi_p = 0.0;
    
    // Dot product µ * r
    double mu2_r = 0.0;
    // Dyadic product Q : R
    double Q2__R = 0.0;

    if (p2._rank > 0) {
        mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
        if (p2._rank > 1) {
            Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                               +   p2.Qzz*rzz);
        }
    }
    
    phi_p = p2.Q00*B0 + mu2_r*B1 + Q2__R*B2;    
    p1.PhiP += phi_p;    
    return phi_p;
}


inline double EwdInteractor::PhiU12_ERFC_At_By(APolarSite &p1, APolarSite &p2) {
    // ATTENTION Potentials are currently not damped
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    
    double phi_u = 0.0;
    
    // Dot product µ(ind) * r
    double u_mu2_r = 0.0;
    
    u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
    
    phi_u = u_mu2_r*B1;
    p1.PhiU += phi_u;
    return phi_u;
}


inline double EwdInteractor::PhiPU12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllCls();
    
    double phi_p = 0.0;
    double phi_u = 0.0;
    
    if (R1 < 1e-2) {
        phi_p += 2.   *a1*EWD::rSqrtPi * (p2.Q00);
        phi_u += 0.0;
    }
    else {
        // Dot product µ * r
        double mu2_r = 0.0;
        // Dot product µ(ind) * r
        double u_mu2_r = 0.0;
        // Dyadic product Q : R
        double Q2__R = 0.0;

        u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);
        if (p2._rank > 0) {
            mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
            if (p2._rank > 1) {
                Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                    +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                   +   p2.Qzz*rzz);
            }
        }

        phi_p = p2.Q00*C0 + mu2_r*C1 + Q2__R*C2;
        phi_u = u_mu2_r*C1;    
    }
    
    // Note the (-): This is a compensation term.
    p1.PhiP -= phi_p;
    p1.PhiU -= phi_u;    
    return phi_p + phi_u;
}


inline double EwdInteractor::PhiP12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllCls();
    
    double phi_p = 0.0;
    
    if (R1 < 1e-2) {
        phi_p += 2.   *a1*EWD::rSqrtPi * (p2.Q00);
    }
    else {
        // Dot product µ * r
        double mu2_r = 0.0;
        // Dyadic product Q : R
        double Q2__R = 0.0;

        if (p2._rank > 0) {
            mu2_r = (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz);
            if (p2._rank > 1) {
                Q2__R = (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                    +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                                   +   p2.Qzz*rzz);
            }
        }

        phi_p = p2.Q00*C0 + mu2_r*C1 + Q2__R*C2;
    }
    
    // Note the (-): This is a compensation term.
    p1.PhiP -= phi_p;
    return phi_p;
}


inline double EwdInteractor::PhiU12_ERF_At_By(APolarSite &p1, APolarSite &p2) {
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    ApplyBiasPolar(p1, p2);
    UpdateAllCls();
    
    double phi_u = 0.0;
    
    if (R1 < 1e-2) {
        phi_u += 0.0;
    }
    else {
        // Dot product µ(ind) * r
        double u_mu2_r = 0.0;

        u_mu2_r = (p2.U1x*rx + p2.U1y*ry + p2.U1z*rz);

        phi_u = u_mu2_r*C1;
    }
    
    // Note the (-): This is a compensation term.
    p1.PhiU -= phi_u;    
    return phi_u;
}


// =============================== REAL SPACE =============================== //
//                                  ENERGIES                                  //

inline EWD::triple<double> EwdInteractor::U12_ERFC(APolarSite &p1, 
    APolarSite &p2) {
    
    ApplyBiasPolar(p1, p2);
    UpdateAllBls();
    UpdateAllGls(p1, p2);
    
    if (R1 < 1e-1) {
        std::cout << std::endl << "small small " << p1.getPos() << " == " << p2.getPos() << flush;
    }       
    
    double pp = 
        ppG0*B0 + ppG1*B1 + ppG2*B2 + ppG3*B3 + ppG4*B4;
    double pu = (ta1*tu3 < 40) ?
        puG1*lp3*B1 + puG2*lp5*B2 + puG3*lp7*B3
      : puG1*B1 + puG2*B2 + puG3*B3;
    double uu = (ta1*tu3 < 40) ?
        uuG1*l3*B1 + uuG2*l5*B2
      : uuG1*B1 + uuG2*B2;
    
    return EWD::triple<double>(pp,pu,uu);
}


inline EWD::triple<double> EwdInteractor::U12_ERF(APolarSite &p1, 
    APolarSite &p2) {
    // NOTE This is a compensation term, it will later be *subtracted* from
    //      the total energy
    
    ApplyBias(p1, p2);
    
    double pp = 0.0;
    double pu = 0.0;
    double uu = 0.0;
    
    if (R1 < 1e-2) {
        //cout << std::endl << "small small " << p1.getPos() << " == " << p2.getPos() << flush;
        pp += 2.   *a1*EWD::rSqrtPi * (p1.Q00*p2.Q00);
        uu += 4./3.*a3*EWD::rSqrtPi * (p1.U1x*p2.U1x + p1.U1y*p2.U1y + p1.U1z*p2.U1z);
        if (p1._rank > 0 && p2._rank > 0) {
            pp += 4./3.*a3*EWD::rSqrtPi * (p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z);
            if (p1._rank > 1 && p2._rank > 1) {
                pp += 16./5.*a5*EWD::rSqrtPi * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
                                                         +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
                                                                           +   p1.Qzz*p2.Qzz);
            }
        }
        if (p1._rank > 0) {
            pu += 4./3.*a3*EWD::rSqrtPi * (p1.Q1x*p2.U1x + p1.Q1y*p2.U1y + p1.Q1z*p2.U1z);
        }
        if (p2._rank > 0) {
            pu += 4./3.*a3*EWD::rSqrtPi * (p1.U1x*p2.Q1x + p1.U1y*p2.Q1y + p1.U1z*p2.Q1z);
        }
        
//        u12 =  2.   *a1*EWD::rSqrtPi * (p1.Q00*p2.Q00)
//            +  4./3.*a3*EWD::rSqrtPi * (p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z)
//            + 16./5.*a5*EWD::rSqrtPi * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
//                                                 +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
//                                                                   +   p1.Qzz*p2.Qzz);
    }
    
    else {
        UpdateAllCls();
        UpdateAllGls(p1, p2);
        pp = ppG0*C0 + ppG1*C1 + ppG2*C2 + ppG3*C3 + ppG4*C4;
        // Damping functions needed here? No.
        pu = puG1*C1 + puG2*C2 + puG3*C3;
        uu = uuG1*C1 + uuG2*C2;
    }
    
    return EWD::triple<double>(pp,pu,uu);
}


//inline EWD::triple<double> EwdInteractor::U12_XYSlab(APolarSite& p1, 
//    APolarSite& p2) {
//    // NOTE : w/o prefactor -2PI/V
//    
//    ApplyBias(p1, p2);    
//    double pp = 0.0;
//    double pu = 0.0;
//    double uu = 0.0;
//    
//    // 1 q <> 2 q
//    pp += p1.Q00*p2.Q00*rz*rz;
//    pu +=  2 * p2.Q00 * p1.U1z * rz;
//    pu += -2 * p1.Q00 * p2.U1z * rz;
//    uu += -2 * p1.U1z * p2.U1z;
//    // 1 d <> 2 q
//    if (p1._rank > 0) {
//        pp += 2 * p2.Q00 * p1.Q1z * rz;
//        pu += -2 * p1.Q1z * p2.U1z;
//        // 1 d <> 2 d
//        if (p2._rank > 0) {
//            pp += -2 * p1.Q1z * p2.Q1z;
//        }
//        // 1 Q <> 2 q
//        if (p1._rank > 1) {
//            pp += 2 * p2.Q00 * p1.Qzz;
//        }
//    }    
//    // 2 d <> 1 q
//    if (p2._rank > 0) {
//        pp += -2 * p1.Q00 * p2.Q1z * rz;
//        pu += -2 * p1.U1z * p2.Q1z;
//        // 2 d <> 1 d
//        // ... covered above.        
//        // 2 Q <> 1 q
//        if (p2._rank > 1) {
//            pp += 2 * p1.Q00 * p2.Qzz;
//        }
//    }
//    
//    return EWD::triple<double>(pp,pu,uu);
//}


inline double EwdInteractor::gG0(APolarSite& p1, APolarSite& p2) {
    return p1.Q00 * p2.Q00;
}


inline double EwdInteractor::gG1(APolarSite& p1, APolarSite& p2) {
    return p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z   
       +   p1.Q00 * (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz)   
       -   p2.Q00 * (p1.Q1x*rx + p1.Q1y*ry + p1.Q1z*rz);
}


inline double EwdInteractor::gG2(APolarSite& p1, APolarSite& p2) {
    return 2 * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
                              +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
                                                +   p1.Qzz*p2.Qzz)
            
       +   p1.Q00 * (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                                +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                               +   p2.Qzz*rzz)
       +   p2.Q00 * (p1.Qxx*rxx + 2*p1.Qxy*rxy + 2*p1.Qxz*rxz
                                +   p1.Qyy*ryy + 2*p1.Qyz*ryz
                                               +   p1.Qzz*rzz)
            
       -   (p1.Q1x*rx + p1.Q1y*ry + p1.Q1z*rz)
          *(p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz)
            
       -   2 * (p1.Qxx*p2.Q1x*rx + p1.Qxy*p2.Q1x*ry + p1.Qxz*p2.Q1x*rz
              + p1.Qxy*p2.Q1y*rx + p1.Qyy*p2.Q1y*ry + p1.Qyz*p2.Q1y*rz
              + p1.Qxz*p2.Q1z*rx + p1.Qyz*p2.Q1z*ry + p1.Qzz*p2.Q1z*rz)
       +   2 * (p2.Qxx*p1.Q1x*rx + p2.Qxy*p1.Q1x*ry + p2.Qxz*p1.Q1x*rz
              + p2.Qxy*p1.Q1y*rx + p2.Qyy*p1.Q1y*ry + p2.Qyz*p1.Q1y*rz
              + p2.Qxz*p1.Q1z*rx + p2.Qyz*p1.Q1z*ry + p2.Qzz*p1.Q1z*rz);
}


inline double EwdInteractor::gG3(APolarSite& p1, APolarSite& p2) {
    return -4 * ( (p1.Qxx*rx + p1.Qxy*ry + p1.Qxz*rz)*(p2.Qxx*rx + p2.Qxy*ry + p2.Qxz*rz)
                + (p1.Qxy*rx + p1.Qyy*ry + p1.Qyz*rz)*(p2.Qxy*rx + p2.Qyy*ry + p2.Qyz*rz)
                + (p1.Qxz*rx + p1.Qyz*ry + p1.Qzz*rz)+(p2.Qxz*rx + p2.Qyz*ry + p2.Qzz*rz) )
            
      -    (p1.Q1x*rx + p1.Q1y*ry + p1.Q1z*rz)
         * (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                      +   p2.Qzz*rzz)
      +    (p2.Q1x*rx + p2.Q1y*ry + p2.Q1z*rz)
         * (p1.Qxx*rxx + 2*p1.Qxy*rxy + 2*p1.Qxz*rxz
                       +   p1.Qyy*ryy + 2*p1.Qyz*ryz
                                      +   p1.Qzz*rzz);
}


inline double EwdInteractor::gG4(APolarSite& p1, APolarSite& p2) {
    return (p1.Qxx*rxx + 2*p1.Qxy*rxy + 2*p1.Qxz*rxz
                       +   p1.Qyy*ryy + 2*p1.Qyz*ryz
                                      +   p1.Qzz*rzz)
         * (p2.Qxx*rxx + 2*p2.Qxy*rxy + 2*p2.Qxz*rxz
                       +   p2.Qyy*ryy + 2*p2.Qyz*ryz
                                      +   p2.Qzz*rzz);
}


// ============================ RECIPROCAL SPACE ============================ //
//                               BIAS FUCTIONS                                //

inline void EwdInteractor::ApplyBiasK(const vec &k) {
    
    k12 = k;
    
    kx = k12.getX(); 
    ky = k12.getY(); 
    kz = k12.getZ();
    
    kxx = kx*kx; kxy = kx*ky; kxz = kx*kz;
    kyy = ky*ky; kyz = ky*kz;
    kzz = kz*kz;
    
    K = votca::tools::abs(k12);
    AK = 4*M_PI*exp(-K*K/(4*a2))/(K*K);
    
    return;
}


inline double EwdInteractor::Ark2Expk2(const vec &k) {
    // 4*M_PI*exp(-K*K/(4*a2))/(K*K)
    ApplyBiasK(k);
    return AK;
}


inline void EwdInteractor::PUApplyBiasK(APolarSite &p) {
    // Permanent & induced moments
    u_dk = p.U1x*kx + p.U1y*ky + p.U1z*kz;
    
    if (p._rank > 0) {
        dk = p.Q1x*kx + p.Q1y*ky + p.Q1z*kz;
        if (p._rank > 1) {
            Qk = p.Qxx*kxx + 2*p.Qxy*kxy + 2*p.Qxz*kxz
                           +   p.Qyy*kyy + 2*p.Qyz*kyz
                                         +   p.Qzz*kzz;
        }
        else Qk = 0.0;
    }        
    else {
        dk = 0.0;
        Qk = 0.0;
    }
    
    kr = kx*p.getPos().getX()
       + ky*p.getPos().getY()
       + kz*p.getPos().getZ();
    coskr = cos(kr);
    sinkr = sin(kr);
    
    re_s   = (p.Q00 - Qk) * coskr   -   dk * sinkr;
    u_re_s =                        - u_dk * sinkr;
    im_s   = (p.Q00 - Qk) * sinkr   +   dk * coskr;
    u_im_s =                        + u_dk * coskr;
    return;
}


inline void EwdInteractor::PApplyBiasK(APolarSite &p) {
    // Permanent moments
    if (p._rank > 0) {
        dk = p.Q1x*kx + p.Q1y*ky + p.Q1z*kz;
        if (p._rank > 1) {
            Qk = p.Qxx*kxx + 2*p.Qxy*kxy + 2*p.Qxz*kxz
                           +   p.Qyy*kyy + 2*p.Qyz*kyz
                                         +   p.Qzz*kzz;
        }
        else Qk = 0.0;
    }        
    else {
        dk = 0.0;
        Qk = 0.0;
    }
    
    kr = kx*p.getPos().getX()
       + ky*p.getPos().getY()
       + kz*p.getPos().getZ();
    coskr = cos(kr);
    sinkr = sin(kr);
    
    re_s   = (p.Q00 - Qk) * coskr   -   dk * sinkr;
    im_s   = (p.Q00 - Qk) * sinkr   +   dk * coskr;
    return;
}


inline void EwdInteractor::UApplyBiasK(APolarSite &p) {
    // Induced moments
    u_dk = p.U1x*kx + p.U1y*ky + p.U1z*kz;
    
    kr = kx*p.getPos().getX()
       + ky*p.getPos().getY()
       + kz*p.getPos().getZ();
    coskr = cos(kr);
    sinkr = sin(kr);
    
    u_re_s =                        - u_dk * sinkr;
    u_im_s =                        + u_dk * coskr;
    return;
}

}}

#endif // VOTCA_XTP_EWDINTERACTOR_H


