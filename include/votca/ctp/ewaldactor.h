#ifndef VOTCA_CTP_EWDINTERACTOR_H
#define	VOTCA_CTP_EWDINTERACTOR_H

#include <cmath>
#include <votca/tools/vec.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/polartop.h>


namespace votca { namespace ctp {

// UNITS IN INPUT FILES
// ... ... Positions as required by format
// ... ... Multipole moment of rank k in e(a0)**k
// ... ... Dipole polarizability in A³ (Angstrom cubed)

// UNITS USED INTERNALLY
// ... ... Positions in nm
// ... ... Multipole moment of rank k in e(nm)**k
// ... ... Dipole polarizability in nm³

// CONVERSION FACTORS
// ... ... Electric field (N/C) = Field (int)  * 1/4PiEps0(SI) * e * 1e+18
// ... ... Energy (eV)          = Energy (int) * 1/4PiEps0(SI) * e * 1e+09
// ... ... Potential(V)         = Pot. (int)   * 1/4PiEps0(SI) * e * 1e+09
    

class EwdInteractor
{
public:

    EwdInteractor(double alpha) {
        a1 = alpha;
        a2 = a1*a1;
        a3 = a1*a2;
        a4 = a1*a3;
        a5 = a1*a4;
        a6 = a1*a5;
        a7 = a1*a6;
    };
    
    EwdInteractor() {};
   ~EwdInteractor() {};
   
    static const double int2eV  = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    static const double rSqrtPi = 0.564189583547756279280349644977832;
    
    struct cmplx
    {
        cmplx(double re, double im) : _re(re), _im(im) { ; }
        double _re;
        double _im;
    };
    
    // ============================= REAL SPACE ============================= //
    
    inline void ApplyBias(APolarSite &p1, APolarSite &p2);
    
    // Make sure to set R1, R2, ... and rR1, rR2, ... before using {gB0, ...}
    inline void UpdateAllBls();
    inline double gB0() { return erfc(a1*R1)*rR1; }
    inline double gB1() { return rR2*(   gB0()  +  2*a1*rSqrtPi * exp(-a2*R2)); }
    inline double gB2() { return rR2*( 3*gB1()  +  4*a3*rSqrtPi * exp(-a2*R2)); }
    inline double gB3() { return rR2*( 5*gB2()  +  8*a5*rSqrtPi * exp(-a2*R2)); }
    inline double gB4() { return rR2*( 7*gB3()  + 16*a7*rSqrtPi * exp(-a2*R2)); }
    
    // Make sure to set R1, R2, ... and rR1, rR2, ... before using {gC0, ...}
    inline void UpdateAllCls();
    inline double gC0() { return erf(a1*R1)*rR1; }
    inline double gC1() { return rR2*(   gC0()  -  2*a1*rSqrtPi * exp(-a2*R2)); }
    inline double gC2() { return rR2*( 3*gC1()  -  4*a3*rSqrtPi * exp(-a2*R2)); }
    inline double gC3() { return rR2*( 5*gC2()  -  8*a5*rSqrtPi * exp(-a2*R2)); }
    inline double gC4() { return rR2*( 7*gC3()  - 16*a7*rSqrtPi * exp(-a2*R2)); }
    
    // Make sure to set R1, R2, ... and rxx, rxy, ... before using {gG0, ...}
    inline void UpdateAllGls(APolarSite &p1, APolarSite &p2);
    inline double gG0(APolarSite &p1, APolarSite &p2);
    inline double gG1(APolarSite &p1, APolarSite &p2);
    inline double gG2(APolarSite &p1, APolarSite &p2);
    inline double gG3(APolarSite &p1, APolarSite &p2);
    inline double gG4(APolarSite &p1, APolarSite &p2);
        
    // Real-space term contribution P1 <> P2
    inline double U12_ERFC(APolarSite &p1, APolarSite &p2);
    
    // Reciprocal-space double-counting correction term P1 <> P2
    inline double U12_ERF(APolarSite &p1, APolarSite &p2);
    
    // Reciprocal-space K=0 shape correction term P1 <> P2
    inline double U12_XYSlab(APolarSite &p1, APolarSite &p2);
    
    
    
    // ========================== RECIPROCAL SPACE ========================== //
    
    inline void ApplyBiasK(const vec &k);
    inline void ApplyBiasK(APolarSite &p);
    inline cmplx AS1S2(const vec &k, vector<PolarSeg*> &s1, vector<PolarSeg*> &s2);
    inline cmplx S1S2(const vec &k, vector<PolarSeg*> &s1, vector<PolarSeg*> &s2);
    inline double Ark2Expk2(const vec &k);
    
    
private:
    
    // Ewald sharpness parameter powers
    double a1, a2, a3, a4, a5, a6, a7;
    
    // ============================= REAL SPACE ============================= //
    
    // Connection vector (1) <- (2), i.e. r12 = r1 - r2 == rab = ra - rb;
    vec r12;
    
    // Vector components rx = r12x, ...
    double rx, ry, rz;
    
    // Matrix product rxx = r12x*r12x, rxy = r12x*r12y, ...
    double rxx, rxy, rxz, ryy, ryz, rzz;
    
    // Real-space distance powers
    double R1, R2, R3, R4, R5;
    
    // Real-space inverse distance powers
    double rR1, rR2, rR3, rR4, rR5;
    
    // {G}, {H} function values
    double G0, G1, G2, G3, G4;
    double H0, H1, H2, H3, H4;
    
    // {Bl}, {Cl} function values
    double rSqrtPiExp;
    double B0, B1, B2, B3, B4;
    double C0, C1, C2, C3, C4;
    
    
    // ========================== RECIPROCAL SPACE ========================== //
    
    // k-space vector
    vec k12;
    double K;
    double AK;    
    
    // Vector components kx, ...
    double kx, ky, kz;
    
    // Matrix product kxx = kx*kx, kxy = kx*ky, ...
    double kxx, kxy, kxz, kyy, kyz, kzz;
    
    // cos(k*r), sin(k*r), µ * k, Q : K
    double dk;
    double Qk;
    double kr;
    double coskr;
    double sinkr;    
    double re_s, im_s;
};


// ============================ RECIPROCAL SPACE ============================ //

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


inline void EwdInteractor::ApplyBiasK(APolarSite &p) {
    
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
    
    re_s = (p.Q00 - Qk) * coskr   -   dk * sinkr;
    im_s = (p.Q00 - Qk) * sinkr   +   dk * coskr;
    return;
}


inline EwdInteractor::cmplx EwdInteractor::AS1S2(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2) {
    // NOTE : w/o 1/V
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // Structure amplitude S1
    double re_S1 = 0.0;
    double im_S1 = 0.0;    
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            ApplyBiasK(*(*pit));            
            re_S1 += re_s;
            im_S1 += im_s; // NOTE THE (+)
        }
    }    
    
    // Structure amplitude S2
    double re_S2 = 0.0;
    double im_S2 = 0.0;    
    for (sit = s2.begin(); sit < s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            ApplyBiasK(*(*pit));
            re_S2 += re_s;
            im_S2 -= im_s; // NOTE THE (-)            
        }
    }
    
    double re_AS1S2 = AK * (re_S1*re_S2 - im_S1*im_S2);
    double im_AS1S2 = AK * (re_S1*im_S2 + im_S1*re_S2);
    
    return cmplx(re_AS1S2, im_AS1S2);
}


inline double EwdInteractor::Ark2Expk2(const vec &k) {
    ApplyBiasK(k);
    return AK;
}


inline EwdInteractor::cmplx EwdInteractor::S1S2(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2) {
    // NOTE : w/o 1/V
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // Structure amplitude S1
    double re_S1 = 0.0;
    double im_S1 = 0.0;    
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            ApplyBiasK(*(*pit));            
            re_S1 += re_s;
            im_S1 += im_s; // NOTE THE (+)
        }
    }    
    
    // Structure amplitude S2
    double re_S2 = 0.0;
    double im_S2 = 0.0;    
    for (sit = s2.begin(); sit < s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            ApplyBiasK(*(*pit));
            re_S2 += re_s;
            im_S2 -= im_s; // NOTE THE (-)            
        }
    }
    
    double re_S1S2 = (re_S1*re_S2 - im_S1*im_S2);
    double im_S1S2 = (re_S1*im_S2 + im_S1*re_S2);
    
    return cmplx(re_S1S2, im_S1S2);
}


// =============================== REAL SPACE =============================== //

inline double EwdInteractor::U12_ERFC(APolarSite &p1, APolarSite &p2) {
    
    ApplyBias(p1, p2);
    UpdateAllBls();
    UpdateAllGls(p1, p2);
    
    if (R1 < 1e-1) {
        cout << endl << "small small " << p1.getPos() << " == " << p2.getPos() << flush;
    }       
    
    return G0*B0 + G1*B1 + G2*B2 + G3*B3 + G4*B4;    
}


inline double EwdInteractor::U12_ERF(APolarSite &p1, APolarSite &p2) {
    
    ApplyBias(p1, p2);
    
    double u12 = 0.0;
    
    if (R1 < 1e-2) {
        //cout << endl << "small small " << p1.getPos() << " == " << p2.getPos() << flush;
        u12 += 2.   *a1*rSqrtPi * (p1.Q00*p2.Q00);
        if (p1._rank > 0 && p2._rank > 0) {
            u12 += 4./3.*a3*rSqrtPi * (p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z);
            if (p1._rank > 1 && p2._rank > 1) {
                u12 += 16./5.*a5*rSqrtPi * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
                                                          +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
                                                                            +   p1.Qzz*p2.Qzz);
            }
        }
        
//        u12 =  2.   *a1*rSqrtPi * (p1.Q00*p2.Q00)
//            +  4./3.*a3*rSqrtPi * (p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z)
//            + 16./5.*a5*rSqrtPi * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
//                                                 +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
//                                                                   +   p1.Qzz*p2.Qzz);
    }
    
    else {
        UpdateAllCls();
        UpdateAllGls(p1, p2);
        u12 = G0*C0 + G1*C1 + G2*C2 + G3*C3 + G4*C4;
    }
    
    return u12;
}


inline double EwdInteractor::U12_XYSlab(APolarSite& p1, APolarSite& p2) {
    // NOTE : w/o prefactor -2PI/V
    
    ApplyBias(p1, p2);    
    double u12 = 0.0;
    
    // 1 q <> 2 q
    u12 += p1.Q00*p2.Q00*rz*rz;    
    // 1 d <> 2 q
    if (p1._rank > 0) {
        u12 += 2 * p2.Q00 * p1.Q1z * rz;        
        // 1 d <> 2 d
        if (p2._rank > 0) {
            u12 += -2 * p1.Q1z * p2.Q1z;
        }
        // 1 Q <> 2 q
        if (p1._rank > 1) {
            u12 += 2 * p2.Q00 * p1.Qzz;
        }
    }    
    // 2 d <> 1 q
    if (p2._rank > 0) {
        u12 += -2 * p1.Q00 * p2.Q1z * rz;        
        // 2 d <> 1 d
        // ... covered above.        
        // 2 Q <> 1 q
        if (p2._rank > 1) {
            u12 += 2 * p1.Q00 * p2.Qzz;
        }
    }
    
    return u12;
}


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


inline void EwdInteractor::UpdateAllBls() {
    
    rSqrtPiExp = rSqrtPi * exp(-a2*R2);
    
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
    
    double rSqrtPiExp = rSqrtPi * exp(-a2*R2);
    
//    C0 = erf(a1*R1)*rR1;    
//    C1 = rR2*(   C0  -  2*a1*rSqrtPiExp);
//    C2 = rR2*( 3*C1  -  4*a3*rSqrtPiExp);
//    C3 = rR2*( 5*C2  -  8*a5*rSqrtPiExp);
//    C4 = rR2*( 7*C3  - 16*a7*rSqrtPiExp);

    C0 = gC0();
    C1 = gC1();
    C2 = gC2();
    C3 = gC3();
    C4 = gC4();
    
    return;
}


inline void EwdInteractor::UpdateAllGls(APolarSite& p1, APolarSite& p2) {
    
    G0 = G1 = G2 = G3 = G4 = 0.0;
    
    // Dot product µ * r
    double mu1_r = 0.0;
    double mu2_r = 0.0;
    
    // Dot product Q * r
    double Q1_rx, Q1_ry, Q1_rz = 0.0;
    double Q2_rx, Q2_ry, Q2_rz = 0.0;
    
    // Dyadic product Q : R
    double Q1__R = 0.0;
    double Q2__R = 0.0;
    
    if (p1._rank > 0) {        
        mu1_r = (p1.Q1x*rx + p1.Q1y*ry + p1.Q1z*rz);        
        if (p1._rank > 1) {
            Q1__R = (p1.Qxx*rxx + 2*p1.Qxy*rxy + 2*p1.Qxz*rxz
                                +   p1.Qyy*ryy + 2*p1.Qyz*ryz
                                               +   p1.Qzz*rzz);
            Q1_rx = p1.Qxx*rx + p1.Qxy*ry + p1.Qxz*rz;
            Q1_ry = p1.Qxy*rx + p1.Qyy*ry + p1.Qyz*rz;
            Q1_rz = p1.Qxz*rx * p1.Qyz*ry + p1.Qzz*rz;
        }
    }
    
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
    
    // 1 - charge, 2 -charge
    G0 = p1.Q00 * p2.Q00;
    
    // 1 - dipole, 2 - charge
    if (p1._rank > 0) {
        G1 += - p2.Q00 * mu1_r;
        
        // 1 - dipole, 2 - dipole
        if (p2._rank > 0) {
            G1 += p1.Q1x*p2.Q1x + p1.Q1y*p2.Q1y + p1.Q1z*p2.Q1z;
            G2 += - mu1_r * mu2_r;
        }
        
        // 1 - quadrupole, 2 -charge
        if (p1._rank > 1) {
            G2 += p2.Q00 * Q1__R;
            
            // 1 - quadrupole, 2 - dipole
            if (p2._rank > 0) {
                G2 += - 2 * (p1.Qxx*p2.Q1x*rx + p1.Qxy*p2.Q1x*ry + p1.Qxz*p2.Q1x*rz
                           + p1.Qxy*p2.Q1y*rx + p1.Qyy*p2.Q1y*ry + p1.Qyz*p2.Q1y*rz
                           + p1.Qxz*p2.Q1z*rx + p1.Qyz*p2.Q1z*ry + p1.Qzz*p2.Q1z*rz);
                G3 += mu2_r * Q1__R;
                
                // 1 - quadrupole, 2 - quadrupole
                if (p2._rank > 1) {
                    G2 += 2 * (p1.Qxx*p2.Qxx + 2*p1.Qxy*p2.Qxy + 2*p1.Qxz*p2.Qxz
                                             +   p1.Qyy*p2.Qyy + 2*p1.Qyz*p2.Qyz
                                                               +   p1.Qzz*p2.Qzz);
                    G3 += -4 * (Q1_rx*Q2_rx + Q1_ry*Q2_ry + Q1_rz*Q2_rz);
                    G4 += Q1__R * Q2__R;
                }
            }            
        }        
    }
    
    // 2 - dipole, 1 - charge
    if (p2._rank > 0) {
        G1 += + p1.Q00 * mu2_r;        
        
        // 2 - dipole, 1 - dipole
        // ... covered above.
        
        // 2 - quadrupole, 1 - charge
        if (p2._rank > 1) {
            G2 += p1.Q00 * Q2__R;
            
            // 2 - quadrupole, 1 - dipole
            if (p1._rank > 0) {
                G2 +=   2 * (p2.Qxx*p1.Q1x*rx + p2.Qxy*p1.Q1x*ry + p2.Qxz*p1.Q1x*rz
                           + p2.Qxy*p1.Q1y*rx + p2.Qyy*p1.Q1y*ry + p2.Qyz*p1.Q1y*rz
                           + p2.Qxz*p1.Q1z*rx + p2.Qyz*p1.Q1z*ry + p2.Qzz*p1.Q1z*rz);
                G3 += - mu1_r * Q2__R;
                
                // 2 - quadrupole, 2 - quadrupole
                // ... covered above.
            }
        }
    }
    
    return;    
}


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


}}

#endif

