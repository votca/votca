/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __XINTERACTOR__H
#define	__XINTERACTOR__H

#include <math.h>
#include <votca/tools/vec.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/apolarsite.h>

// TODO Sharpness parameter should be initialised in constructor
// ...  (currently hard-coded, 0.390)


namespace votca { namespace ctp {

// +++++++++++++++++++++++++++ //
// Multipole Interaction Class //
// +++++++++++++++++++++++++++ //

class XInteractor
{
public:

    XInteractor(Topology *top, double aDamp) : _top(top),  a(aDamp) {};
    XInteractor() :                            _top(NULL), a(0.390) {};
   ~XInteractor() {};

    // UNITS IN INPUT FILES
    // ... Always use atomic units
    // ... ... Positions in a0 (bohr)
    // ... ... Multipole moment of rank k in e(a0)**k
    // ... ... Dipole polarizability in A³ (Angstrom cubed)

    // UNITS USED INTERNALLY
    // ... Use nm instead of a0 and A
    // ... ... Positions in nm
    // ... ... Multipole moment of rank k in e(nm)**k
    // ... ... Dipole polarizability in nm³

    // CONVERSION FACTORS
    // ... Electric field (N/C) = Electric field (int)
    //                          * 1/4PiEps0(SI) * e * 1e+18
    // ... Energy (eV) = Energy (int) * 1/4PiEps0(SI) * e * 1e+09
    //
    // ... Potential (V) = Potential(int) * 1/4PiEps0(SI) * e * 1e+0.9


    inline double   EnergyInter(APolarSite &pol1, APolarSite &pol2);
    inline double   EnergyInterESP(APolarSite &pol1, APolarSite &pol2);
    inline double   EnergyIntra(APolarSite &pol1, APolarSite &pol2);
    inline void     FieldPerm(APolarSite &pol1, APolarSite &pol2);
    inline vec      FieldPermESF(vec r, APolarSite &pol);
    inline void     FieldIndu(APolarSite &pol1, APolarSite &pol2);
    inline void     FieldInduAlpha(APolarSite &pol1, APolarSite &pol2);

    inline double   PotentialPerm(vec r, APolarSite &pol);

    inline double   E_f(APolarSite &pol1, APolarSite &pol2);
    inline double   E_m(APolarSite &pol1, APolarSite &pol2);
    
    inline double   E_QQ_ERFC(APolarSite &pol1, APolarSite &pol2, double &ew_alpha);
    inline double   E_QQ_ERF(APolarSite &pol1, APolarSite &pol2, double &ew_alpha);
    inline double   E_QQ_K0(APolarSite &pol1, APolarSite &pol2, double &ew_alpha);
    inline double   E_QQ_KK(APolarSite &pol1, APolarSite &pol2, double &ew_alpha, vec &k);
    inline double   E_Q0_DQ(APolarSite &pol1, APolarSite &pol2);

    void            ResetEnergy() { EP = EU_INTER = EU_INTRA = 0.0;
                                    EPP = EPU = EUU = 0.0; }
    double          &getEP() { return EP; }
    double          &getEU_INTER() { return EU_INTER; }
    double          &getEU_INTRA() { return EU_INTRA; }

    double          &getEPP() { return EPP; }
    double          &getEPU() { return EPU; }
    double          &getEUU() { return EUU; }
    
//    inline void     FieldIndu12_21(vector<PolarSite*> &, vector<PolarSite*> &);
//    inline void     FieldIndu12___(vector<PolarSite*> &, vector<PolarSite*> &);
//    inline void     FieldIndu___21(vector<PolarSite*> &, vector<PolarSite*> &);
    
    inline void     BiasIndu(APolarSite &pol1, APolarSite &pol2);
    inline void     BiasStat(APolarSite &pol1, APolarSite &pol2);
    inline void     RevBias();
    

private:

    double EP;       //   <- Interaction permanent multipoles (inter-site)
    double EU_INTRA; //   <- Interaction induction multipoles (intra-site)
    double EU_INTER; //   <- Interaction induction multipoles (inter-site)

    double EPP;
    double EPU;
    double EUU;


    vec    e12;     //  |
    double u3;      //  |-> NOTE: Only needed when using Thole model
    double a;       //  |         (do not forget to init. though...)

    double R;       //  |
    double R2;      //  |
    double R3;      //  |-> NOTE: reciprocal, i.e. e.g. R3 = 1/(R*R*R)
    double R4;      //  |
    double R5;      //  |

    double rax, ray, raz;
    double rbx, rby, rbz;
    double cxx, cxy, cxz;
    double cyx, cyy, cyz;
    double czx, czy, czz;

    inline double lambda3() { return 1 - exp( -a*u3); }
    inline double lambda5() { return 1 - (1 + a*u3) * exp( -a*u3); }
    inline double lambda7() { return 1 - (1 + a*u3 + 0.6*a*a*u3*u3) * exp( -a*u3); }
    inline double lambda9() { return 1 - (1 + a*u3 + (18*a*a*u3*u3 + 9*a*a*a*u3*u3*u3)/35) * exp( -a*u3); }

    inline double T00_00() { return R; }

    inline double T1x_00() { return R2 * rax; }
    inline double T1y_00() { return R2 * ray; }
    inline double T1z_00() { return R2 * raz; }
    inline double T00_1x() { return R2 * rbx; }
    inline double T00_1y() { return R2 * rby; }
    inline double T00_1z() { return R2 * rbz; }

    inline double TU1x_00() { return lambda3() * R2 * rax; }
    inline double TU1y_00() { return lambda3() * R2 * ray; }
    inline double TU1z_00() { return lambda3() * R2 * raz; }
    inline double TU00_1x() { return lambda3() * R2 * rbx; }
    inline double TU00_1y() { return lambda3() * R2 * rby; }
    inline double TU00_1z() { return lambda3() * R2 * rbz; }

    inline double T20_00()  { return R3 * 0.5 * (3 * raz*raz - 1); }
    inline double T21c_00() { return R3 * sqrt(3) * rax * raz; }
    inline double T21s_00() { return R3 * sqrt(3) * ray * raz; }
    inline double T22c_00() { return R3 * 0.5 * sqrt(3) * (rax*rax - ray*ray); }
    inline double T22s_00() { return R3 * sqrt(3) * rax*ray; }
    inline double T00_20()  { return R3 * 0.5 * (3 * rbz*rbz - 1); }
    inline double T00_21c() { return R3 * sqrt(3) * rbx * rbz; }
    inline double T00_21s() { return R3 * sqrt(3) * rby * rbz; }
    inline double T00_22c() { return R3 * 0.5 * sqrt(3) * (rbx*rbx - rby*rby); }
    inline double T00_22s() { return R3 * sqrt(3) * rbx*rby; }

    inline double T1x_1x() { return R3 * (3 * rax*rbx + cxx); }
    inline double T1x_1y() { return R3 * (3 * rax*rby + cxy); }
    inline double T1x_1z() { return R3 * (3 * rax*rbz + cxz); }
    inline double T1y_1x() { return R3 * (3 * ray*rbx + cyx); }
    inline double T1y_1y() { return R3 * (3 * ray*rby + cyy); }
    inline double T1y_1z() { return R3 * (3 * ray*rbz + cyz); }
    inline double T1z_1x() { return R3 * (3 * raz*rbx + czx); }
    inline double T1z_1y() { return R3 * (3 * raz*rby + czy); }
    inline double T1z_1z() { return R3 * (3 * raz*rbz + czz); }

    inline double TU1x_1x() { return R3 * (lambda5()*3*rax*rbx + lambda3()*cxx); }
    inline double TU1x_1y() { return R3 * (lambda5()*3*rax*rby + lambda3()*cxy); }
    inline double TU1x_1z() { return R3 * (lambda5()*3*rax*rbz + lambda3()*cxz); }
    inline double TU1y_1x() { return R3 * (lambda5()*3*ray*rbx + lambda3()*cyx); }
    inline double TU1y_1y() { return R3 * (lambda5()*3*ray*rby + lambda3()*cyy); }
    inline double TU1y_1z() { return R3 * (lambda5()*3*ray*rbz + lambda3()*cyz); }
    inline double TU1z_1x() { return R3 * (lambda5()*3*raz*rbx + lambda3()*czx); }
    inline double TU1z_1y() { return R3 * (lambda5()*3*raz*rby + lambda3()*czy); }
    inline double TU1z_1z() { return R3 * (lambda5()*3*raz*rbz + lambda3()*czz); }

    inline double T20_1x()  { return R4 * 0.5 * (15*raz*raz*rbx + 6*raz*czx - 3*rbx); }
    inline double T20_1y()  { return R4 * 0.5 * (15*raz*raz*rby + 6*raz*czy - 3*rby); }
    inline double T20_1z()  { return R4 * 0.5 * (15*raz*raz*rbz + 6*raz*czz - 3*rbz); }
    inline double T21c_1x() { return R4 * sqrt(3) * (rax*czx + cxx*raz + 5*rax*raz*rbx); }
    inline double T21c_1y() { return R4 * sqrt(3) * (rax*czy + cxy*raz + 5*rax*raz*rby); }
    inline double T21c_1z() { return R4 * sqrt(3) * (rax*czz + cxz*raz + 5*rax*raz*rbz); }
    inline double T21s_1x() { return R4 * sqrt(3) * (ray*czx + cyx*raz + 5*ray*raz*rbx); }
    inline double T21s_1y() { return R4 * sqrt(3) * (ray*czy + cyy*raz + 5*ray*raz*rby); }
    inline double T21s_1z() { return R4 * sqrt(3) * (ray*czz + cyz*raz + 5*ray*raz*rbz); }
    inline double T22c_1x() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rbx + 2*rax*cxx - 2*ray*cyx); }
    inline double T22c_1y() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rby + 2*rax*cxy - 2*ray*cyy); }
    inline double T22c_1z() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rbz + 2*rax*cxz - 2*ray*cyz); }
    inline double T22s_1x() { return R4 * sqrt(3) * ( 5*rax*ray*rbx + rax*cyx + ray*cxx ); }
    inline double T22s_1y() { return R4 * sqrt(3) * ( 5*rax*ray*rby + rax*cyy + ray*cxy ); }
    inline double T22s_1z() { return R4 * sqrt(3) * ( 5*rax*ray*rbz + rax*cyz + ray*cxz ); }

    inline double T1x_20()  { return R4 * 0.5 * (15*rbz*rbz*rax + 6*rbz*cxz - 3*rax); }
    inline double T1y_20()  { return R4 * 0.5 * (15*rbz*rbz*ray + 6*rbz*cyz - 3*ray); }
    inline double T1z_20()  { return R4 * 0.5 * (15*rbz*rbz*raz + 6*rbz*czz - 3*raz); }
    inline double T1x_21c() { return R4 * sqrt(3) * (rbx*cxz + cxx*rbz + 5*rbx*rbz*rax); }
    inline double T1y_21c() { return R4 * sqrt(3) * (rbx*cyz + cyx*rbz + 5*rbx*rbz*ray); }
    inline double T1z_21c() { return R4 * sqrt(3) * (rbx*czz + czx*rbz + 5*rbx*rbz*raz); }
    inline double T1x_21s() { return R4 * sqrt(3) * (rby*cxz + cxy*rbz + 5*rby*rbz*rax); }
    inline double T1y_21s() { return R4 * sqrt(3) * (rby*cyz + cyy*rbz + 5*rby*rbz*ray); }
    inline double T1z_21s() { return R4 * sqrt(3) * (rby*czz + czy*rbz + 5*rby*rbz*raz); }
    inline double T1x_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*rax + 2*rbx*cxx - 2*rby*cxy); }
    inline double T1y_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*ray + 2*rbx*cyx - 2*rby*cyy); }
    inline double T1z_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*raz + 2*rbx*czx - 2*rby*czy); }
    inline double T1x_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*rax + rbx*cxy + rby*cxx ); }
    inline double T1y_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*ray + rbx*cyy + rby*cyx ); }
    inline double T1z_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*raz + rbx*czy + rby*czx ); }

    inline double TU20_1x()  { return R4 * 0.5 * (lambda7()*15*raz*raz*rbx + lambda5()*(6*raz*czx - 3*rbx)); }
    inline double TU20_1y()  { return R4 * 0.5 * (lambda7()*15*raz*raz*rby + lambda5()*(6*raz*czy - 3*rby)); }
    inline double TU20_1z()  { return R4 * 0.5 * (lambda7()*15*raz*raz*rbz + lambda5()*(6*raz*czz - 3*rbz)); }
    inline double TU21c_1x() { return R4 * sqrt(3) * (lambda5()*(rax*czx + cxx*raz) + lambda7()*5*rax*raz*rbx); }
    inline double TU21c_1y() { return R4 * sqrt(3) * (lambda5()*(rax*czy + cxy*raz) + lambda7()*5*rax*raz*rby); }
    inline double TU21c_1z() { return R4 * sqrt(3) * (lambda5()*(rax*czz + cxz*raz) + lambda7()*5*rax*raz*rbz); }
    inline double TU21s_1x() { return R4 * sqrt(3) * (lambda5()*(ray*czx + cyx*raz) + lambda7()*5*ray*raz*rbx); }
    inline double TU21s_1y() { return R4 * sqrt(3) * (lambda5()*(ray*czy + cyy*raz) + lambda7()*5*ray*raz*rby); }
    inline double TU21s_1z() { return R4 * sqrt(3) * (lambda5()*(ray*czz + cyz*raz) + lambda7()*5*ray*raz*rbz); }
    inline double TU22c_1x() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rax*rax-ray*ray)*rbx + lambda5()*(2*rax*cxx - 2*ray*cyx)); }
    inline double TU22c_1y() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rax*rax-ray*ray)*rby + lambda5()*(2*rax*cxy - 2*ray*cyy)); }
    inline double TU22c_1z() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rax*rax-ray*ray)*rbz + lambda5()*(2*rax*cxz - 2*ray*cyz)); }
    inline double TU22s_1x() { return R4 * sqrt(3) * (lambda7()*5*rax*ray*rbx + lambda5()*(rax*cyx + ray*cxx) ); }
    inline double TU22s_1y() { return R4 * sqrt(3) * (lambda7()*5*rax*ray*rby + lambda5()*(rax*cyy + ray*cxy) ); }
    inline double TU22s_1z() { return R4 * sqrt(3) * (lambda7()*5*rax*ray*rbz + lambda5()*(rax*cyz + ray*cxz) ); }

    inline double TU1x_20()  { return R4 * 0.5 * (lambda7()*15*rbz*rbz*rax + lambda5()*(6*rbz*cxz - 3*rax)); }
    inline double TU1y_20()  { return R4 * 0.5 * (lambda7()*15*rbz*rbz*ray + lambda5()*(6*rbz*cyz - 3*ray)); }
    inline double TU1z_20()  { return R4 * 0.5 * (lambda7()*15*rbz*rbz*raz + lambda5()*(6*rbz*czz - 3*raz)); }
    inline double TU1x_21c() { return R4 * sqrt(3) * (lambda5()*(rbx*cxz + cxx*rbz) + lambda7()*5*rbx*rbz*rax); }
    inline double TU1y_21c() { return R4 * sqrt(3) * (lambda5()*(rbx*cyz + cyx*rbz) + lambda7()*5*rbx*rbz*ray); }
    inline double TU1z_21c() { return R4 * sqrt(3) * (lambda5()*(rbx*czz + czx*rbz) + lambda7()*5*rbx*rbz*raz); }
    inline double TU1x_21s() { return R4 * sqrt(3) * (lambda5()*(rby*cxz + cxy*rbz) + lambda7()*5*rby*rbz*rax); }
    inline double TU1y_21s() { return R4 * sqrt(3) * (lambda5()*(rby*cyz + cyy*rbz) + lambda7()*5*rby*rbz*ray); }
    inline double TU1z_21s() { return R4 * sqrt(3) * (lambda5()*(rby*czz + czy*rbz) + lambda7()*5*rby*rbz*raz); }
    inline double TU1x_22c() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rbx*rbx-rby*rby)*rax + lambda5()*(2*rbx*cxx - 2*rby*cxy)); }
    inline double TU1y_22c() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rbx*rbx-rby*rby)*ray + lambda5()*(2*rbx*cyx - 2*rby*cyy)); }
    inline double TU1z_22c() { return R4 * 0.5 * sqrt(3) * (lambda7()*5*(rbx*rbx-rby*rby)*raz + lambda5()*(2*rbx*czx - 2*rby*czy)); }
    inline double TU1x_22s() { return R4 * sqrt(3) * (lambda7()*5*rbx*rby*rax + lambda5()*(rbx*cxy + rby*cxx) ); }
    inline double TU1y_22s() { return R4 * sqrt(3) * (lambda7()*5*rbx*rby*ray + lambda5()*(rbx*cyy + rby*cyx) ); }
    inline double TU1z_22s() { return R4 * sqrt(3) * (lambda7()*5*rbx*rby*raz + lambda5()*(rbx*czy + rby*czx) ); }

    inline double T20_20()   { return R5 * 0.75 * (35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz + 20*raz*rbz*czz + 2*czz*czz + 1); }
    inline double T20_21c()  { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rbx*rbz - 5*rbx*rbz + 10*raz*rbx*czz + 10*raz*rbz*czx + 2*czx*czz); }
    inline double T20_21s()  { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rby*rbz - 5*rby*rbz + 10*raz*rby*czz + 10*raz*rbz*czy + 2*czy*czz); }
    inline double T20_22c()  { return R5 * 0.25 * sqrt(3) * (35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby - 5*rbx*rbx + 5*rby*rby + 20*raz*rbx*czx - 20*raz*rby*czy + 2*czx*czx - 2*czy*czy); }
    inline double T20_22s()  { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rbx*rby - 5*rbx*rby + 10*raz*rbx*czy + 10*raz*rby*czx + 2*czx*czy); }
    inline double T21c_21c() { return R5 * (35*rax*raz*rbx*rbz + 5*rax*rbx*czz + 5*rax*rbz*czx + 5*raz*rbx*cxz + 5*raz*rbz*cxx + cxx*czz + cxz*czx); }
    inline double T21c_21s() { return R5 * (35*rax*raz*rby*rbz + 5*rax*rby*czz + 5*rax*rbz*czy + 5*raz*rby*cxz + 5*raz*rbz*cxy + cxy*czz + cxz*czy); }
    inline double T21c_22c() { return R5 * 0.5 * (35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby + 10*rax*rbx*czx - 10*rax*rby*czy + 10*raz*rbx*cxx - 10*raz*rby*cxy + 2*cxx*czx - 2*cxy*czy); }
    inline double T21c_22s() { return R5 * (35*rax*raz*rbx*rby + 5*rax*rbx*czy + 5*rax*rby*czx + 5*raz*rbx*cxy + 5*raz*rby*cxx + cxx*czy + cxy*czx); }
    inline double T21s_21s() { return R5 * (35*ray*raz*rby*rbz + 5*ray*rby*czz + 5*ray*rbz*czy + 5*raz*rby*cyz + 5*raz*rbz*cyy + cyy*czz + cyz*czy); }
    inline double T21s_22c() { return R5 * 0.5 * (35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby + 10*ray*rbx*czx - 10*ray*rby*czy + 10*raz*rbx*cyx - 10*raz*rby*cyy + 2*cyx*czx - 2*cyy*czy); }
    inline double T21s_22s() { return R5 * (35*ray*raz*rbx*rby + 5*ray*rbx*czy + 5*ray*rby*czx + 5*raz*rbx*cyy + 5*raz*rby*cyx + cyx*czy + cyy*czx); }
    inline double T22c_22c() { return R5 * 0.25 * (35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby + 20*rax*rbx*cxx - 20*rax*rby*cxy - 20*ray*rbx*cyx + 20*ray*rby*cyy + 2*cxx*cxx - 2*cxy*cxy - 2*cyx*cyx + 2*cyy*cyy); }
    inline double T22c_22s() { return R5 * 0.5 * (35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby + 10*rax*rbx*cxy + 10*rax*rby*cxx - 10*ray*rbx*cyy - 10*ray*rby*cyx + 2*cxx*cxy - 2*cyx*cyy); }
    inline double T22s_22s() { return R5 * (35*rax*ray*rbx*rby + 5*rax*rbx*cyy + 5*rax*rby*cyx + 5*ray*rbx*cxy + 5*ray*rby*cxx + cxx*cyy + cxy*cyx); }

    inline double T21c_20()  { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*rax*raz - 5*rax*raz + 10*rbz*rax*czz + 10*rbz*raz*cxz + 2*cxz*czz); }
    inline double T21s_20()  { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*ray*raz - 5*ray*raz + 10*rbz*ray*czz + 10*rbz*raz*cyz + 2*cyz*czz); }
    inline double T22c_20()  { return R5 * 0.25 * sqrt(3) * (35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray - 5*rax*rax + 5*ray*ray + 20*rbz*rax*cxz - 20*rbz*ray*cyz + 2*cxz*cxz - 2*cyz*cyz); }
    inline double T22s_20()  { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*rax*ray - 5*rax*ray + 10*rbz*rax*cyz + 10*rbz*ray*cxz + 2*cxz*cyz); }
    inline double T21s_21c() { return R5 * (35*rbx*rbz*ray*raz + 5*rbx*ray*czz + 5*rbx*raz*cyz + 5*rbz*ray*czx + 5*rbz*raz*cyx + cyx*czz + czx*cyz); }
    inline double T22c_21c() { return R5 * 0.5 * (35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray + 10*rbx*rax*cxz - 10*rbx*ray*cyz + 10*rbz*rax*cxx - 10*rbz*ray*cyx + 2*cxx*cxz - 2*cyx*cyz); }
    inline double T22s_21c() { return R5 * (35*rbx*rbz*rax*ray + 5*rbx*rax*cyz + 5*rbx*ray*cxz + 5*rbz*rax*cyx + 5*rbz*ray*cxx + cxx*cyz + cyx*cxz); }
    inline double T22c_21s() { return R5 * 0.5 * (35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray + 10*rby*rax*cxz - 10*rby*ray*cyz + 10*rbz*rax*cxy - 10*rbz*ray*cyy + 2*cxy*cxz - 2*cyy*cyz); }
    inline double T22s_21s() { return R5 * (35*rby*rbz*rax*ray + 5*rby*rax*cyz + 5*rby*ray*cxz + 5*rbz*rax*cyy + 5*rbz*ray*cxy + cxy*cyz + cyy*cxz); }
    inline double T22s_22c() { return R5 * 0.5 * (35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray + 10*rbx*rax*cyx + 10*rbx*ray*cxx - 10*rby*rax*cyy - 10*rby*ray*cxy + 2*cxx*cyx - 2*cxy*cyy); }

    Topology        *_top;

};

    
    
// ========================================================================== //
//                         XINTERACTOR MEMBER FUNCTIONS                       //
// ========================================================================== //


/**
 * Used in ESP calculator (initialize stage of XQMP)
 */
inline double XInteractor::PotentialPerm(vec r, APolarSite &pol) {

    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol.getPos() - r;
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    


    rbx = - e12.getX();
    rby = - e12.getX();
    rbz = - e12.getX();

    double phi00 = 0.0;
    
        phi00 += T00_00() * pol.Q00;
    
    if (pol._rank > 0) {        
        phi00 += T00_1x() * pol.Q1x;
        phi00 += T00_1y() * pol.Q1y;
        phi00 += T00_1z() * pol.Q1z;        
    }
        
    if (pol._rank > 1) {
        phi00 += T00_20()  * pol.Q20;
        phi00 += T00_21c() * pol.Q21c;
        phi00 += T00_21s() * pol.Q21s;
        phi00 += T00_22c() * pol.Q22c;
        phi00 += T00_22s() * pol.Q22s;
    }

    return phi00;
}

/**
 * Used in ESF calculator (initialize stage of XQMP)
 */
inline vec XInteractor::FieldPermESF(vec r, APolarSite &pol) {
    
    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol.getPos() - r;
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;

    double Fx = 0.0;
    double Fy = 0.0;
    double Fz = 0.0;

    // Field generated by rank-0 m'pole
        Fx += T1x_00() * pol.Q00;
        Fy += T1y_00() * pol.Q00;
        Fz += T1z_00() * pol.Q00;

    // Field generated by rank-1 m'pole
    if (pol._rank > 0) {
        Fx += T1x_1x() * pol.Q1x;
        Fx += T1x_1y() * pol.Q1y;
        Fx += T1x_1z() * pol.Q1z;

        Fy += T1y_1x() * pol.Q1x;
        Fy += T1y_1y() * pol.Q1y;
        Fy += T1y_1z() * pol.Q1z;

        Fz += T1z_1x() * pol.Q1x;
        Fz += T1z_1y() * pol.Q1y;
        Fz += T1z_1z() * pol.Q1z;
    }

    // Field generated by rank-2 m'pole
    if (pol._rank > 1) {
        Fx += T1x_20()  * pol.Q20;
        Fx += T1x_21c() * pol.Q21c;
        Fx += T1x_21s() * pol.Q21s;
        Fx += T1x_22c() * pol.Q22c;
        Fx += T1x_22s() * pol.Q22s;

        Fy += T1y_20()  * pol.Q20;
        Fy += T1y_21c() * pol.Q21c;
        Fy += T1y_21s() * pol.Q21s;
        Fy += T1y_22c() * pol.Q22c;
        Fy += T1y_22s() * pol.Q22s;

        Fz += T1z_20()  * pol.Q20;
        Fz += T1z_21c() * pol.Q21c;
        Fz += T1z_21s() * pol.Q21s;
        Fz += T1z_22c() * pol.Q22c;
        Fz += T1z_22s() * pol.Q22s;
    }

    return vec(Fx, Fy, Fz);
}

/**
 * Used in molecular-polarizability calculator (initialize stage)
 */
inline void XInteractor::FieldInduAlpha(APolarSite &pol1, APolarSite &pol2) {
    
    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));

//        rax =   pol1._locX * e12;
//        ray =   pol1._locY * e12;
//        raz =   pol1._locZ * e12;
//        rbx = - pol2._locX * e12;
//        rby = - pol2._locY * e12;
//        rbz = - pol2._locZ * e12;

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

//        cxx = pol1._locX * pol2._locX;
//        cxy = pol1._locX * pol2._locY;
//        cxz = pol1._locX * pol2._locZ;
//        cyx = pol1._locY * pol2._locX;
//        cyy = pol1._locY * pol2._locY;
//        cyz = pol1._locY * pol2._locZ;
//        czx = pol1._locZ * pol2._locX;
//        czy = pol1._locZ * pol2._locY;
//        czz = pol1._locZ * pol2._locZ;

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;

    // Fields generated by rank-1 induced m'poles

    if (a*u3 < 40.0) {
        pol1.FUx += TU1x_1x() * pol2.U1x;
        pol1.FUx += TU1x_1y() * pol2.U1y;
        pol1.FUx += TU1x_1z() * pol2.U1z;
        pol1.FUy += TU1y_1x() * pol2.U1x;
        pol1.FUy += TU1y_1y() * pol2.U1y;
        pol1.FUy += TU1y_1z() * pol2.U1z;
        pol1.FUz += TU1z_1x() * pol2.U1x;
        pol1.FUz += TU1z_1y() * pol2.U1y;
        pol1.FUz += TU1z_1z() * pol2.U1z;

        pol2.FUx += TU1x_1x() * pol1.U1x;
        pol2.FUx += TU1y_1x() * pol1.U1y;
        pol2.FUx += TU1z_1x() * pol1.U1z;
        pol2.FUy += TU1x_1y() * pol1.U1x;
        pol2.FUy += TU1y_1y() * pol1.U1y;
        pol2.FUy += TU1z_1y() * pol1.U1z;
        pol2.FUz += TU1x_1z() * pol1.U1x;
        pol2.FUz += TU1y_1z() * pol1.U1y;
        pol2.FUz += TU1z_1z() * pol1.U1z;
    }
    else {
        pol1.FUx += T1x_1x() * pol2.U1x;
        pol1.FUx += T1x_1y() * pol2.U1y;
        pol1.FUx += T1x_1z() * pol2.U1z;
        pol1.FUy += T1y_1x() * pol2.U1x;
        pol1.FUy += T1y_1y() * pol2.U1y;
        pol1.FUy += T1y_1z() * pol2.U1z;
        pol1.FUz += T1z_1x() * pol2.U1x;
        pol1.FUz += T1z_1y() * pol2.U1y;
        pol1.FUz += T1z_1z() * pol2.U1z;

        pol2.FUx += T1x_1x() * pol1.U1x;
        pol2.FUx += T1y_1x() * pol1.U1y;
        pol2.FUx += T1z_1x() * pol1.U1z;
        pol2.FUy += T1x_1y() * pol1.U1x;
        pol2.FUy += T1y_1y() * pol1.U1y;
        pol2.FUy += T1z_1y() * pol1.U1z;
        pol2.FUz += T1x_1z() * pol1.U1x;
        pol2.FUz += T1y_1z() * pol1.U1y;
        pol2.FUz += T1z_1z() * pol1.U1z;
    }
}

/**
 * Used in self-consistent field calculation (evaluation stage)
 */
inline void XInteractor::FieldIndu(APolarSite &pol1, APolarSite &pol2) {

//    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
//    //          This implies that induced = - alpha * field
//    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
//    //e12  = pol2.getPos() - pol1.getPos();
//    R    = 1/abs(e12);
//    R2   = R*R;
//    R3   = R2*R;
//    R4   = R3*R;
//    R5   = R4*R;
//    e12 *= R;
//
//    // Thole damping init.
//    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));
//
////        rax =   pol1._locX * e12;
////        ray =   pol1._locY * e12;
////        raz =   pol1._locZ * e12;
////        rbx = - pol2._locX * e12;
////        rby = - pol2._locY * e12;
////        rbz = - pol2._locZ * e12;
//
//        rax = e12.getX();
//        ray = e12.getY();
//        raz = e12.getZ();
//        rbx = - rax;
//        rby = - ray;
//        rbz = - raz;
//
////        cxx = pol1._locX * pol2._locX;
////        cxy = pol1._locX * pol2._locY;
////        cxz = pol1._locX * pol2._locZ;
////        cyx = pol1._locY * pol2._locX;
////        cyy = pol1._locY * pol2._locY;
////        cyz = pol1._locY * pol2._locZ;
////        czx = pol1._locZ * pol2._locX;
////        czy = pol1._locZ * pol2._locY;
////        czz = pol1._locZ * pol2._locZ;
//
//        cxx = 1;
//        cxy = 0;
//        cxz = 0;
//        cyx = 0;
//        cyy = 1;
//        cyz = 0;
//        czx = 0;
//        czy = 0;
//        czz = 1;

    // Fields generated by rank-1 induced m'poles

    if (a*u3 < 40.0) {
        pol1.FUx += TU1x_1x() * pol2.U1x;
        pol1.FUx += TU1x_1y() * pol2.U1y;
        pol1.FUx += TU1x_1z() * pol2.U1z;
        pol1.FUy += TU1y_1x() * pol2.U1x;
        pol1.FUy += TU1y_1y() * pol2.U1y;
        pol1.FUy += TU1y_1z() * pol2.U1z;
        pol1.FUz += TU1z_1x() * pol2.U1x;
        pol1.FUz += TU1z_1y() * pol2.U1y;
        pol1.FUz += TU1z_1z() * pol2.U1z;

        pol2.FUx += TU1x_1x() * pol1.U1x;
        pol2.FUx += TU1y_1x() * pol1.U1y;
        pol2.FUx += TU1z_1x() * pol1.U1z;
        pol2.FUy += TU1x_1y() * pol1.U1x;
        pol2.FUy += TU1y_1y() * pol1.U1y;
        pol2.FUy += TU1z_1y() * pol1.U1z;
        pol2.FUz += TU1x_1z() * pol1.U1x;
        pol2.FUz += TU1y_1z() * pol1.U1y;
        pol2.FUz += TU1z_1z() * pol1.U1z;
    }
    else {
        pol1.FUx += T1x_1x() * pol2.U1x;
        pol1.FUx += T1x_1y() * pol2.U1y;
        pol1.FUx += T1x_1z() * pol2.U1z;
        pol1.FUy += T1y_1x() * pol2.U1x;
        pol1.FUy += T1y_1y() * pol2.U1y;
        pol1.FUy += T1y_1z() * pol2.U1z;
        pol1.FUz += T1z_1x() * pol2.U1x;
        pol1.FUz += T1z_1y() * pol2.U1y;
        pol1.FUz += T1z_1z() * pol2.U1z;

        pol2.FUx += T1x_1x() * pol1.U1x;
        pol2.FUx += T1y_1x() * pol1.U1y;
        pol2.FUx += T1z_1x() * pol1.U1z;
        pol2.FUy += T1x_1y() * pol1.U1x;
        pol2.FUy += T1y_1y() * pol1.U1y;
        pol2.FUy += T1z_1y() * pol1.U1z;
        pol2.FUz += T1x_1z() * pol1.U1x;
        pol2.FUz += T1y_1z() * pol1.U1y;
        pol2.FUz += T1z_1z() * pol1.U1z;  
    }
}

/**
 * Used in self-consistent field calculation (evaluation stage)
 */
inline void XInteractor::FieldPerm(APolarSite &pol1, APolarSite &pol2) {
    
//    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
//    //          This implies that induced = - alpha * field
//    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
//    //e12  = pol2.getPos() - pol1.getPos();
//    R    = 1/abs(e12);
//    R2   = R*R;
//    R3   = R2*R;
//    R4   = R3*R;
//    R5   = R4*R;
//    e12 *= R;
//
////        rax =   pol1._locX * e12;
////        ray =   pol1._locY * e12;
////        raz =   pol1._locZ * e12;
////        rbx = - pol2._locX * e12;
////        rby = - pol2._locY * e12;
////        rbz = - pol2._locZ * e12;
//
//        rax = e12.getX();
//        ray = e12.getY();
//        raz = e12.getZ();
//        rbx = - rax;
//        rby = - ray;
//        rbz = - raz;
//
//    if (pol1._rank > 0 || pol2._rank > 0) {
////        cxx = pol1._locX * pol2._locX;
////        cxy = pol1._locX * pol2._locY;
////        cxz = pol1._locX * pol2._locZ;
////        cyx = pol1._locY * pol2._locX;
////        cyy = pol1._locY * pol2._locY;
////        cyz = pol1._locY * pol2._locZ;
////        czx = pol1._locZ * pol2._locX;
////        czy = pol1._locZ * pol2._locY;
////        czz = pol1._locZ * pol2._locZ;
//
//        cxx = 1;
//        cxy = 0;
//        cxz = 0;
//        cyx = 0;
//        cyy = 1;
//        cyz = 0;
//        czx = 0;
//        czy = 0;
//        czz = 1;
//    }

    // Fields generated by rank-0 m'poles
        pol1.FPx += T1x_00() * pol2.Q00;
        pol1.FPy += T1y_00() * pol2.Q00;
        pol1.FPz += T1z_00() * pol2.Q00;

        pol2.FPx += T00_1x() * pol1.Q00;
        pol2.FPy += T00_1y() * pol1.Q00;
        pol2.FPz += T00_1z() * pol1.Q00;

    // Fields generated by rank-1 m'poles
    if (pol2._rank > 0) {
        pol1.FPx += T1x_1x() * pol2.Q1x;
        pol1.FPx += T1x_1y() * pol2.Q1y;
        pol1.FPx += T1x_1z() * pol2.Q1z;
        pol1.FPy += T1y_1x() * pol2.Q1x;
        pol1.FPy += T1y_1y() * pol2.Q1y;
        pol1.FPy += T1y_1z() * pol2.Q1z;
        pol1.FPz += T1z_1x() * pol2.Q1x;
        pol1.FPz += T1z_1y() * pol2.Q1y;
        pol1.FPz += T1z_1z() * pol2.Q1z;
    }
    if (pol1._rank > 0) {
        pol2.FPx += T1x_1x() * pol1.Q1x;
        pol2.FPx += T1y_1x() * pol1.Q1y;
        pol2.FPx += T1z_1x() * pol1.Q1z;
        pol2.FPy += T1x_1y() * pol1.Q1x;
        pol2.FPy += T1y_1y() * pol1.Q1y;
        pol2.FPy += T1z_1y() * pol1.Q1z;
        pol2.FPz += T1x_1z() * pol1.Q1x;
        pol2.FPz += T1y_1z() * pol1.Q1y;
        pol2.FPz += T1z_1z() * pol1.Q1z;
    }

    // Fields generated by rank-2 m'poles
    if (pol2._rank > 1) {
        pol1.FPx += T1x_20()  * pol2.Q20;
        pol1.FPx += T1x_21c() * pol2.Q21c;
        pol1.FPx += T1x_21s() * pol2.Q21s;
        pol1.FPx += T1x_22c() * pol2.Q22c;
        pol1.FPx += T1x_22s() * pol2.Q22s;

        pol1.FPy += T1y_20()  * pol2.Q20;
        pol1.FPy += T1y_21c() * pol2.Q21c;
        pol1.FPy += T1y_21s() * pol2.Q21s;
        pol1.FPy += T1y_22c() * pol2.Q22c;
        pol1.FPy += T1y_22s() * pol2.Q22s;

        pol1.FPz += T1z_20()  * pol2.Q20;
        pol1.FPz += T1z_21c() * pol2.Q21c;
        pol1.FPz += T1z_21s() * pol2.Q21s;
        pol1.FPz += T1z_22c() * pol2.Q22c;
        pol1.FPz += T1z_22s() * pol2.Q22s;
    }
    if (pol1._rank > 1) {
        pol2.FPx += T20_1x()  * pol1.Q20;
        pol2.FPx += T21c_1x() * pol1.Q21c;
        pol2.FPx += T21s_1x() * pol1.Q21s;
        pol2.FPx += T22c_1x() * pol1.Q22c;
        pol2.FPx += T22s_1x() * pol1.Q22s;

        pol2.FPy += T20_1y()  * pol1.Q20;
        pol2.FPy += T21c_1y() * pol1.Q21c;
        pol2.FPy += T21s_1y() * pol1.Q21s;
        pol2.FPy += T22c_1y() * pol1.Q22c;
        pol2.FPy += T22s_1y() * pol1.Q22s;

        pol2.FPz += T20_1z()  * pol1.Q20;
        pol2.FPz += T21c_1z() * pol1.Q21c;
        pol2.FPz += T21s_1z() * pol1.Q21s;
        pol2.FPz += T22c_1z() * pol1.Q22c;
        pol2.FPz += T22s_1z() * pol1.Q22s;        
    }
}

/**
 * Used in energy evaluation of converged fields (evaluation stage)
 */
inline double XInteractor::EnergyIntra(APolarSite &pol1, APolarSite &pol2) {    

    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;
    
    // Thole damping init.
    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));

//        rax =   pol1._locX * e12;
//        ray =   pol1._locY * e12;
//        raz =   pol1._locZ * e12;
//        rbx = - pol2._locX * e12;
//        rby = - pol2._locY * e12;
//        rbz = - pol2._locZ * e12;

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

    if (pol1._rank > 0 || pol2._rank > 0) {
//        cxx = pol1._locX * pol2._locX;
//        cxy = pol1._locX * pol2._locY;
//        cxz = pol1._locX * pol2._locZ;
//        cyx = pol1._locY * pol2._locX;
//        cyy = pol1._locY * pol2._locY;
//        cyz = pol1._locY * pol2._locZ;
//        czx = pol1._locZ * pol2._locX;
//        czy = pol1._locZ * pol2._locY;
//        czz = pol1._locZ * pol2._locZ;

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;
    }

    double U = 0.0; // <- Induction energy

    if (a*u3 < 40.0) {
        U += pol1.U1x * TU1x_00() * pol2.Q00;
        U += pol1.U1y * TU1y_00() * pol2.Q00;
        U += pol1.U1z * TU1z_00() * pol2.Q00;

        U += pol1.Q00 * TU00_1x() * pol2.U1x;
        U += pol1.Q00 * TU00_1y() * pol2.U1y;
        U += pol1.Q00 * TU00_1z() * pol2.U1z;
    }
    else {
        U += pol1.U1x * T1x_00() * pol2.Q00;
        U += pol1.U1y * T1y_00() * pol2.Q00;
        U += pol1.U1z * T1z_00() * pol2.Q00;

        U += pol1.Q00 * T00_1x() * pol2.U1x;
        U += pol1.Q00 * T00_1y() * pol2.U1y;
        U += pol1.Q00 * T00_1z() * pol2.U1z;
    }

    if (pol1._rank > 0) {
        if (a*u3 < 40.0) {
            U += pol1.Q1x * TU1x_1x() * pol2.U1x;
            U += pol1.Q1x * TU1x_1y() * pol2.U1y;
            U += pol1.Q1x * TU1x_1z() * pol2.U1z;
            U += pol1.Q1y * TU1y_1x() * pol2.U1x;
            U += pol1.Q1y * TU1y_1y() * pol2.U1y;
            U += pol1.Q1y * TU1y_1z() * pol2.U1z;
            U += pol1.Q1z * TU1z_1x() * pol2.U1x;
            U += pol1.Q1z * TU1z_1y() * pol2.U1y;
            U += pol1.Q1z * TU1z_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q1x * T1x_1x() * pol2.U1x;
            U += pol1.Q1x * T1x_1y() * pol2.U1y;
            U += pol1.Q1x * T1x_1z() * pol2.U1z;
            U += pol1.Q1y * T1y_1x() * pol2.U1x;
            U += pol1.Q1y * T1y_1y() * pol2.U1y;
            U += pol1.Q1y * T1y_1z() * pol2.U1z;
            U += pol1.Q1z * T1z_1x() * pol2.U1x;
            U += pol1.Q1z * T1z_1y() * pol2.U1y;
            U += pol1.Q1z * T1z_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 0) {
        if (a*u3 < 40.0) {
            U += pol1.U1x * TU1x_1x() * pol2.Q1x;
            U += pol1.U1x * TU1x_1y() * pol2.Q1y;
            U += pol1.U1x * TU1x_1z() * pol2.Q1z;
            U += pol1.U1y * TU1y_1x() * pol2.Q1x;
            U += pol1.U1y * TU1y_1y() * pol2.Q1y;
            U += pol1.U1y * TU1y_1z() * pol2.Q1z;
            U += pol1.U1z * TU1z_1x() * pol2.Q1x;
            U += pol1.U1z * TU1z_1y() * pol2.Q1y;
            U += pol1.U1z * TU1z_1z() * pol2.Q1z;
        }
        else {
            U += pol1.U1x * T1x_1x() * pol2.Q1x;
            U += pol1.U1x * T1x_1y() * pol2.Q1y;
            U += pol1.U1x * T1x_1z() * pol2.Q1z;
            U += pol1.U1y * T1y_1x() * pol2.Q1x;
            U += pol1.U1y * T1y_1y() * pol2.Q1y;
            U += pol1.U1y * T1y_1z() * pol2.Q1z;
            U += pol1.U1z * T1z_1x() * pol2.Q1x;
            U += pol1.U1z * T1z_1y() * pol2.Q1y;
            U += pol1.U1z * T1z_1z() * pol2.Q1z;
        }
    }

    if (pol1._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.Q20  * TU20_1x()  * pol2.U1x;
            U += pol1.Q20  * TU20_1y()  * pol2.U1y;
            U += pol1.Q20  * TU20_1z()  * pol2.U1z;
            U += pol1.Q21c * TU21c_1x() * pol2.U1x;
            U += pol1.Q21c * TU21c_1y() * pol2.U1y;
            U += pol1.Q21c * TU21c_1z() * pol2.U1z;
            U += pol1.Q21s * TU21s_1x() * pol2.U1x;
            U += pol1.Q21s * TU21s_1y() * pol2.U1y;
            U += pol1.Q21s * TU21s_1z() * pol2.U1z;
            U += pol1.Q22c * TU22c_1x() * pol2.U1x;
            U += pol1.Q22c * TU22c_1y() * pol2.U1y;
            U += pol1.Q22c * TU22c_1z() * pol2.U1z;
            U += pol1.Q22s * TU22s_1x() * pol2.U1x;
            U += pol1.Q22s * TU22s_1y() * pol2.U1y;
            U += pol1.Q22s * TU22s_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q20  * T20_1x()  * pol2.U1x;
            U += pol1.Q20  * T20_1y()  * pol2.U1y;
            U += pol1.Q20  * T20_1z()  * pol2.U1z;
            U += pol1.Q21c * T21c_1x() * pol2.U1x;
            U += pol1.Q21c * T21c_1y() * pol2.U1y;
            U += pol1.Q21c * T21c_1z() * pol2.U1z;
            U += pol1.Q21s * T21s_1x() * pol2.U1x;
            U += pol1.Q21s * T21s_1y() * pol2.U1y;
            U += pol1.Q21s * T21s_1z() * pol2.U1z;
            U += pol1.Q22c * T22c_1x() * pol2.U1x;
            U += pol1.Q22c * T22c_1y() * pol2.U1y;
            U += pol1.Q22c * T22c_1z() * pol2.U1z;
            U += pol1.Q22s * T22s_1x() * pol2.U1x;
            U += pol1.Q22s * T22s_1y() * pol2.U1y;
            U += pol1.Q22s * T22s_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.U1x * TU1x_20()  * pol2.Q20;
            U += pol1.U1x * TU1x_21c() * pol2.Q21c;
            U += pol1.U1x * TU1x_21s() * pol2.Q21s;
            U += pol1.U1x * TU1x_22c() * pol2.Q22c;
            U += pol1.U1x * TU1x_22s() * pol2.Q22s;
            U += pol1.U1y * TU1y_20()  * pol2.Q20;
            U += pol1.U1y * TU1y_21c() * pol2.Q21c;
            U += pol1.U1y * TU1y_21s() * pol2.Q21s;
            U += pol1.U1y * TU1y_22c() * pol2.Q22c;
            U += pol1.U1y * TU1y_22s() * pol2.Q22s;
            U += pol1.U1z * TU1z_20()  * pol2.Q20;
            U += pol1.U1z * TU1z_21c() * pol2.Q21c;
            U += pol1.U1z * TU1z_21s() * pol2.Q21s;
            U += pol1.U1z * TU1z_22c() * pol2.Q22c;
            U += pol1.U1z * TU1z_22s() * pol2.Q22s;
        }
        else {
            U += pol1.U1x * T1x_20()  * pol2.Q20;
            U += pol1.U1x * T1x_21c() * pol2.Q21c;
            U += pol1.U1x * T1x_21s() * pol2.Q21s;
            U += pol1.U1x * T1x_22c() * pol2.Q22c;
            U += pol1.U1x * T1x_22s() * pol2.Q22s;
            U += pol1.U1y * T1y_20()  * pol2.Q20;
            U += pol1.U1y * T1y_21c() * pol2.Q21c;
            U += pol1.U1y * T1y_21s() * pol2.Q21s;
            U += pol1.U1y * T1y_22c() * pol2.Q22c;
            U += pol1.U1y * T1y_22s() * pol2.Q22s;
            U += pol1.U1z * T1z_20()  * pol2.Q20;
            U += pol1.U1z * T1z_21c() * pol2.Q21c;
            U += pol1.U1z * T1z_21s() * pol2.Q21s;
            U += pol1.U1z * T1z_22c() * pol2.Q22c;
            U += pol1.U1z * T1z_22s() * pol2.Q22s;
        }
    }
    
    // Take into account work needed to induce multipoles
    U *= 0.5;

    EU_INTRA += U;
    return U;
}

/**
 * Used in energy evaluation of converged fields (evaluation stage)
 */
inline double XInteractor::EnergyInter(APolarSite &pol1, APolarSite &pol2) {

    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));


    //cout << "frag1 " << pol1.getFragment()->getId() << endl;
    //cout << "frag2 " << pol2.getFragment()->getId() << endl;
    //cout << "seg1  " << pol1.getSegment()->getId() << endl;
    //cout << "seg2  " << pol2.getSegment()->getId() << endl;    

    
//        rax =   pol1._locX * e12;
//        ray =   pol1._locY * e12;
//        raz =   pol1._locZ * e12;
//        rbx = - pol2._locX * e12;
//        rby = - pol2._locY * e12;
//        rbz = - pol2._locZ * e12;

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

    if (pol1._rank > 0 || pol2._rank > 0) {
//        cxx = pol1._locX * pol2._locX;
//        cxy = pol1._locX * pol2._locY;
//        cxz = pol1._locX * pol2._locZ;
//        cyx = pol1._locY * pol2._locX;
//        cyy = pol1._locY * pol2._locY;
//        cyz = pol1._locY * pol2._locZ;
//        czx = pol1._locZ * pol2._locX;
//        czy = pol1._locZ * pol2._locY;
//        czz = pol1._locZ * pol2._locZ;

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;
    }

    double E = 0.0; // <- Electrostatic energy
    double U = 0.0; // <- Induction energy

        //cout << "r1  " << pol1.getPos() << endl;
        //cout << "r2  " << pol2.getPos() << endl;
        //cout << "R   " << 1/R << endl;
        //cout << "e12 " << e12 << endl;

        E += pol1.Q00 * T00_00() * pol2.Q00;

        //cout << "E up to q <-> q " << E << endl;

    if (a*u3 < 40) {
        U += pol1.U1x * TU1x_00() * pol2.Q00;
        U += pol1.U1y * TU1y_00() * pol2.Q00;
        U += pol1.U1z * TU1z_00() * pol2.Q00;

        U += pol1.Q00 * TU00_1x() * pol2.U1x;
        U += pol1.Q00 * TU00_1y() * pol2.U1y;
        U += pol1.Q00 * TU00_1z() * pol2.U1z;
    }
    else {
        U += pol1.U1x * T1x_00() * pol2.Q00;
        U += pol1.U1y * T1y_00() * pol2.Q00;
        U += pol1.U1z * T1z_00() * pol2.Q00;

        U += pol1.Q00 * T00_1x() * pol2.U1x;
        U += pol1.Q00 * T00_1y() * pol2.U1y;
        U += pol1.Q00 * T00_1z() * pol2.U1z;
    }



    if (pol1._rank > 0) {
        E += pol1.Q1x * T1x_00() * pol2.Q00;
        //cout << "E1x_00 " << pol1.Q1x * T1x_00() * pol2.Q00 << endl;
        E += pol1.Q1y * T1y_00() * pol2.Q00;
        //cout << "E1y_00 " << pol1.Q1y * T1y_00() * pol2.Q00 << endl;
        E += pol1.Q1z * T1z_00() * pol2.Q00;
        //cout << "E1z_00 " << pol1.Q1z * T1z_00() * pol2.Q00 << endl;
    }

    if (pol2._rank > 0) {
        E += pol1.Q00 * T00_1x() * pol2.Q1x;
        //cout << "E00_1x " << pol1.Q00 * T00_1x() * pol2.Q1x << endl;
        E += pol1.Q00 * T00_1y() * pol2.Q1y;
        //cout << "E00_1y " << pol1.Q00 * T00_1y() * pol2.Q1y << endl;
        E += pol1.Q00 * T00_1z() * pol2.Q1z;
        //cout << "E00_1z " << pol1.Q00 * T00_1z() * pol2.Q1z << endl;
    }
        //cout << "E up to q <-> d " << E << endl;

    if (pol1._rank > 1) {
        E += pol1.Q20  * T20_00()  * pol2.Q00;
        E += pol1.Q21c * T21c_00() * pol2.Q00;
        E += pol1.Q21s * T21s_00() * pol2.Q00;
        E += pol1.Q22c * T22c_00() * pol2.Q00;
        E += pol1.Q22s * T22s_00() * pol2.Q00;
    }

    if (pol2._rank > 1) {
        E += pol1.Q00 * T00_20()  * pol2.Q20;
        E += pol1.Q00 * T00_21c() * pol2.Q21c;
        E += pol1.Q00 * T00_21s() * pol2.Q21s;
        E += pol1.Q00 * T00_22c() * pol2.Q22c;
        E += pol1.Q00 * T00_22s() * pol2.Q22s;
    }
        //cout << "E up to q <-> Q " << E << endl;

    if (pol1._rank > 0 && pol2._rank > 0) {
        E += pol1.Q1x * T1x_1x() * pol2.Q1x;
        //cout << "E1x_1x " << pol1.Q1x * T1x_1x() * pol2.Q1x << endl;
        E += pol1.Q1x * T1x_1y() * pol2.Q1y;
        //cout << "E1x_1y " << pol1.Q1x * T1x_1y() * pol2.Q1y << endl;
        E += pol1.Q1x * T1x_1z() * pol2.Q1z;
        //cout << "E1x_1z " << pol1.Q1x * T1x_1z() * pol2.Q1z << endl;

        E += pol1.Q1y * T1y_1x() * pol2.Q1x;
        //cout << "E1y_1x " << pol1.Q1y * T1y_1x() * pol2.Q1x << endl;
        E += pol1.Q1y * T1y_1y() * pol2.Q1y;
        //cout << "E1y_1y " << pol1.Q1y * T1y_1y() * pol2.Q1y << endl;
        E += pol1.Q1y * T1y_1z() * pol2.Q1z;
        //cout << "E1y_1z " << pol1.Q1y * T1y_1z() * pol2.Q1z << endl;

        E += pol1.Q1z * T1z_1x() * pol2.Q1x;
        //cout << "E1z_1x " << pol1.Q1z * T1z_1x() * pol2.Q1x << endl;
        E += pol1.Q1z * T1z_1y() * pol2.Q1y;
        //cout << "E1z_1y " << pol1.Q1z * T1z_1y() * pol2.Q1y << endl;
        E += pol1.Q1z * T1z_1z() * pol2.Q1z;
        //cout << "E1z_1z " << pol1.Q1z * T1z_1z() * pol2.Q1z << endl;
    }

    if (pol1._rank > 0) {
        if (a*u3 < 40) {
            U += pol1.Q1x * TU1x_1x() * pol2.U1x;
            U += pol1.Q1x * TU1x_1y() * pol2.U1y;
            U += pol1.Q1x * TU1x_1z() * pol2.U1z;
            U += pol1.Q1y * TU1y_1x() * pol2.U1x;
            U += pol1.Q1y * TU1y_1y() * pol2.U1y;
            U += pol1.Q1y * TU1y_1z() * pol2.U1z;
            U += pol1.Q1z * TU1z_1x() * pol2.U1x;
            U += pol1.Q1z * TU1z_1y() * pol2.U1y;
            U += pol1.Q1z * TU1z_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q1x * T1x_1x() * pol2.U1x;
            U += pol1.Q1x * T1x_1y() * pol2.U1y;
            U += pol1.Q1x * T1x_1z() * pol2.U1z;
            U += pol1.Q1y * T1y_1x() * pol2.U1x;
            U += pol1.Q1y * T1y_1y() * pol2.U1y;
            U += pol1.Q1y * T1y_1z() * pol2.U1z;
            U += pol1.Q1z * T1z_1x() * pol2.U1x;
            U += pol1.Q1z * T1z_1y() * pol2.U1y;
            U += pol1.Q1z * T1z_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 0) {
        if (a*u3 < 40) {
            U += pol1.U1x * TU1x_1x() * pol2.Q1x;
            U += pol1.U1x * TU1x_1y() * pol2.Q1y;
            U += pol1.U1x * TU1x_1z() * pol2.Q1z;
            U += pol1.U1y * TU1y_1x() * pol2.Q1x;
            U += pol1.U1y * TU1y_1y() * pol2.Q1y;
            U += pol1.U1y * TU1y_1z() * pol2.Q1z;
            U += pol1.U1z * TU1z_1x() * pol2.Q1x;
            U += pol1.U1z * TU1z_1y() * pol2.Q1y;
            U += pol1.U1z * TU1z_1z() * pol2.Q1z;
        }
        else {
            U += pol1.U1x * T1x_1x() * pol2.Q1x;
            U += pol1.U1x * T1x_1y() * pol2.Q1y;
            U += pol1.U1x * T1x_1z() * pol2.Q1z;
            U += pol1.U1y * T1y_1x() * pol2.Q1x;
            U += pol1.U1y * T1y_1y() * pol2.Q1y;
            U += pol1.U1y * T1y_1z() * pol2.Q1z;
            U += pol1.U1z * T1z_1x() * pol2.Q1x;
            U += pol1.U1z * T1z_1y() * pol2.Q1y;
            U += pol1.U1z * T1z_1z() * pol2.Q1z;
        }
    }    
        //cout << "E up to d <-> d " << E << endl;

    if (pol1._rank > 1 && pol2._rank > 0) {
        E += pol1.Q20 * T20_1x() * pol2.Q1x;
        E += pol1.Q20 * T20_1y() * pol2.Q1y;
        E += pol1.Q20 * T20_1z() * pol2.Q1z;

        E += pol1.Q21c * T21c_1x() * pol2.Q1x;
        E += pol1.Q21c * T21c_1y() * pol2.Q1y;
        E += pol1.Q21c * T21c_1z() * pol2.Q1z;

        E += pol1.Q21s * T21s_1x() * pol2.Q1x;
        E += pol1.Q21s * T21s_1y() * pol2.Q1y;
        E += pol1.Q21s * T21s_1z() * pol2.Q1z;

        E += pol1.Q22c * T22c_1x() * pol2.Q1x;
        E += pol1.Q22c * T22c_1y() * pol2.Q1y;
        E += pol1.Q22c * T22c_1z() * pol2.Q1z;

        E += pol1.Q22s * T22s_1x() * pol2.Q1x;
        E += pol1.Q22s * T22s_1y() * pol2.Q1y;
        E += pol1.Q22s * T22s_1z() * pol2.Q1z;
    }

    if (pol1._rank > 0 && pol2._rank > 1) {
        E += pol1.Q1x * T1x_20() * pol2.Q20;
        E += pol1.Q1y * T1y_20() * pol2.Q20;
        E += pol1.Q1z * T1z_20() * pol2.Q20;

        E += pol1.Q1x * T1x_21c() * pol2.Q21c;
        E += pol1.Q1y * T1y_21c() * pol2.Q21c;
        E += pol1.Q1z * T1z_21c() * pol2.Q21c;

        E += pol1.Q1x * T1x_21s() * pol2.Q21s;
        E += pol1.Q1y * T1y_21s() * pol2.Q21s;
        E += pol1.Q1z * T1z_21s() * pol2.Q21s;

        E += pol1.Q1x * T1x_22c() * pol2.Q22c;
        E += pol1.Q1y * T1y_22c() * pol2.Q22c;
        E += pol1.Q1z * T1z_22c() * pol2.Q22c;

        E += pol1.Q1x * T1x_22s() * pol2.Q22s;
        E += pol1.Q1y * T1y_22s() * pol2.Q22s;
        E += pol1.Q1z * T1z_22s() * pol2.Q22s;
    }

    if (pol1._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.Q20  * TU20_1x()  * pol2.U1x;
            U += pol1.Q20  * TU20_1y()  * pol2.U1y;
            U += pol1.Q20  * TU20_1z()  * pol2.U1z;
            U += pol1.Q21c * TU21c_1x() * pol2.U1x;
            U += pol1.Q21c * TU21c_1y() * pol2.U1y;
            U += pol1.Q21c * TU21c_1z() * pol2.U1z;
            U += pol1.Q21s * TU21s_1x() * pol2.U1x;
            U += pol1.Q21s * TU21s_1y() * pol2.U1y;
            U += pol1.Q21s * TU21s_1z() * pol2.U1z;
            U += pol1.Q22c * TU22c_1x() * pol2.U1x;
            U += pol1.Q22c * TU22c_1y() * pol2.U1y;
            U += pol1.Q22c * TU22c_1z() * pol2.U1z;
            U += pol1.Q22s * TU22s_1x() * pol2.U1x;
            U += pol1.Q22s * TU22s_1y() * pol2.U1y;
            U += pol1.Q22s * TU22s_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q20  * T20_1x()  * pol2.U1x;
            U += pol1.Q20  * T20_1y()  * pol2.U1y;
            U += pol1.Q20  * T20_1z()  * pol2.U1z;
            U += pol1.Q21c * T21c_1x() * pol2.U1x;
            U += pol1.Q21c * T21c_1y() * pol2.U1y;
            U += pol1.Q21c * T21c_1z() * pol2.U1z;
            U += pol1.Q21s * T21s_1x() * pol2.U1x;
            U += pol1.Q21s * T21s_1y() * pol2.U1y;
            U += pol1.Q21s * T21s_1z() * pol2.U1z;
            U += pol1.Q22c * T22c_1x() * pol2.U1x;
            U += pol1.Q22c * T22c_1y() * pol2.U1y;
            U += pol1.Q22c * T22c_1z() * pol2.U1z;
            U += pol1.Q22s * T22s_1x() * pol2.U1x;
            U += pol1.Q22s * T22s_1y() * pol2.U1y;
            U += pol1.Q22s * T22s_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.U1x * TU1x_20()  * pol2.Q20;
            U += pol1.U1x * TU1x_21c() * pol2.Q21c;
            U += pol1.U1x * TU1x_21s() * pol2.Q21s;
            U += pol1.U1x * TU1x_22c() * pol2.Q22c;
            U += pol1.U1x * TU1x_22s() * pol2.Q22s;
            U += pol1.U1y * TU1y_20()  * pol2.Q20;
            U += pol1.U1y * TU1y_21c() * pol2.Q21c;
            U += pol1.U1y * TU1y_21s() * pol2.Q21s;
            U += pol1.U1y * TU1y_22c() * pol2.Q22c;
            U += pol1.U1y * TU1y_22s() * pol2.Q22s;
            U += pol1.U1z * TU1z_20()  * pol2.Q20;
            U += pol1.U1z * TU1z_21c() * pol2.Q21c;
            U += pol1.U1z * TU1z_21s() * pol2.Q21s;
            U += pol1.U1z * TU1z_22c() * pol2.Q22c;
            U += pol1.U1z * TU1z_22s() * pol2.Q22s;
        }
        else {
            U += pol1.U1x * T1x_20()  * pol2.Q20;
            U += pol1.U1x * T1x_21c() * pol2.Q21c;
            U += pol1.U1x * T1x_21s() * pol2.Q21s;
            U += pol1.U1x * T1x_22c() * pol2.Q22c;
            U += pol1.U1x * T1x_22s() * pol2.Q22s;
            U += pol1.U1y * T1y_20()  * pol2.Q20;
            U += pol1.U1y * T1y_21c() * pol2.Q21c;
            U += pol1.U1y * T1y_21s() * pol2.Q21s;
            U += pol1.U1y * T1y_22c() * pol2.Q22c;
            U += pol1.U1y * T1y_22s() * pol2.Q22s;
            U += pol1.U1z * T1z_20()  * pol2.Q20;
            U += pol1.U1z * T1z_21c() * pol2.Q21c;
            U += pol1.U1z * T1z_21s() * pol2.Q21s;
            U += pol1.U1z * T1z_22c() * pol2.Q22c;
            U += pol1.U1z * T1z_22s() * pol2.Q22s;
        }
    }
        //cout << "E up to d <-> Q " << E << endl;

    if (pol1._rank > 1 && pol2._rank > 1) {
        E += pol1.Q20  * T20_20()   * pol2.Q20;
        E += pol1.Q21c * T21c_21c() * pol2.Q21c;
        E += pol1.Q21s * T21s_21s() * pol2.Q21s;
        E += pol1.Q22c * T22c_22c() * pol2.Q22c;
        E += pol1.Q22s * T22s_22s() * pol2.Q22s;


        E += pol1.Q20  * T20_21c() * pol2.Q21c;
        E += pol1.Q20  * T20_21s() * pol2.Q21s;
        E += pol1.Q20  * T20_22c() * pol2.Q22c;
        E += pol1.Q20  * T20_22s() * pol2.Q22s;
        E += pol1.Q21c * T21c_20() * pol2.Q20;
        E += pol1.Q21s * T21s_20() * pol2.Q20;
        E += pol1.Q22c * T22c_20() * pol2.Q20;
        E += pol1.Q22s * T22s_20() * pol2.Q20;


        E += pol1.Q21c * T21c_21s() * pol2.Q21s;
        E += pol1.Q21c * T21c_22c() * pol2.Q22c;
        E += pol1.Q21c * T21c_22s() * pol2.Q22s;
        E += pol1.Q21s * T21s_21c() * pol2.Q21c;
        E += pol1.Q22c * T22c_21c() * pol2.Q21c;
        E += pol1.Q22s * T22s_21c() * pol2.Q21c;


        E += pol1.Q21s * T21s_22c() * pol2.Q22c;
        E += pol1.Q21s * T21s_22s() * pol2.Q22s;
        E += pol1.Q22c * T22c_21s() * pol2.Q21s;
        E += pol1.Q22s * T22s_21s() * pol2.Q21s;

        E += pol1.Q22s * T22s_22c() * pol2.Q22c;
        E += pol1.Q22c * T22c_22s() * pol2.Q22s;
    }
        //cout << "E up to Q <-> Q " << E << endl;
        

    // Take into account work required to induce multipoles
    U *= 0.5;

    EP += E;
    EU_INTER += U;
    return E + U;
}

/**
 * Designed for use in ESP calculator (init. stage). Only for error-checking.
 */
inline double XInteractor::EnergyInterESP(APolarSite &pol1, APolarSite &pol2) {
    
    assert(false);
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

    if (pol1._rank > 0 || pol2._rank > 0) {

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;
    }

    double E = 0.0; // <- Electrostatic energy
    double U = 0.0; // <- Induction energy

        E += pol1.Q00 * T00_00() * pol2.Q00;

    if (a*u3 < 40) {
        U += pol1.U1x * TU1x_00() * pol2.Q00;
        U += pol1.U1y * TU1y_00() * pol2.Q00;
        U += pol1.U1z * TU1z_00() * pol2.Q00;

        U += pol1.Q00 * TU00_1x() * pol2.U1x;
        U += pol1.Q00 * TU00_1y() * pol2.U1y;
        U += pol1.Q00 * TU00_1z() * pol2.U1z;
    }
    else {
        U += pol1.U1x * T1x_00() * pol2.Q00;
        U += pol1.U1y * T1y_00() * pol2.Q00;
        U += pol1.U1z * T1z_00() * pol2.Q00;

        U += pol1.Q00 * T00_1x() * pol2.U1x;
        U += pol1.Q00 * T00_1y() * pol2.U1y;
        U += pol1.Q00 * T00_1z() * pol2.U1z;
    }



    if (pol1._rank > 0) {
        E += pol1.Q1x * T1x_00() * pol2.Q00;
        E += pol1.Q1y * T1y_00() * pol2.Q00;
        E += pol1.Q1z * T1z_00() * pol2.Q00;
    }

    if (pol2._rank > 0) {
        E += pol1.Q00 * T00_1x() * pol2.Q1x;
        E += pol1.Q00 * T00_1y() * pol2.Q1y;
        E += pol1.Q00 * T00_1z() * pol2.Q1z;
    }

    if (pol1._rank > 1) {
        E += pol1.Q20  * T20_00()  * pol2.Q00;
        E += pol1.Q21c * T21c_00() * pol2.Q00;
        E += pol1.Q21s * T21s_00() * pol2.Q00;
        E += pol1.Q22c * T22c_00() * pol2.Q00;
        E += pol1.Q22s * T22s_00() * pol2.Q00;
    }

    if (pol2._rank > 1) {
        E += pol1.Q00 * T00_20()  * pol2.Q20;
        E += pol1.Q00 * T00_21c() * pol2.Q21c;
        E += pol1.Q00 * T00_21s() * pol2.Q21s;
        E += pol1.Q00 * T00_22c() * pol2.Q22c;
        E += pol1.Q00 * T00_22s() * pol2.Q22s;
    }

    if (pol1._rank > 0 && pol2._rank > 0) {
        E += pol1.Q1x * T1x_1x() * pol2.Q1x;
        E += pol1.Q1x * T1x_1y() * pol2.Q1y;
        E += pol1.Q1x * T1x_1z() * pol2.Q1z;

        E += pol1.Q1y * T1y_1x() * pol2.Q1x;
        E += pol1.Q1y * T1y_1y() * pol2.Q1y;
        E += pol1.Q1y * T1y_1z() * pol2.Q1z;

        E += pol1.Q1z * T1z_1x() * pol2.Q1x;
        E += pol1.Q1z * T1z_1y() * pol2.Q1y;
        E += pol1.Q1z * T1z_1z() * pol2.Q1z;
    }

    if (pol1._rank > 0) {
        if (a*u3 < 40) {
            U += pol1.Q1x * TU1x_1x() * pol2.U1x;
            U += pol1.Q1x * TU1x_1y() * pol2.U1y;
            U += pol1.Q1x * TU1x_1z() * pol2.U1z;
            U += pol1.Q1y * TU1y_1x() * pol2.U1x;
            U += pol1.Q1y * TU1y_1y() * pol2.U1y;
            U += pol1.Q1y * TU1y_1z() * pol2.U1z;
            U += pol1.Q1z * TU1z_1x() * pol2.U1x;
            U += pol1.Q1z * TU1z_1y() * pol2.U1y;
            U += pol1.Q1z * TU1z_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q1x * T1x_1x() * pol2.U1x;
            U += pol1.Q1x * T1x_1y() * pol2.U1y;
            U += pol1.Q1x * T1x_1z() * pol2.U1z;
            U += pol1.Q1y * T1y_1x() * pol2.U1x;
            U += pol1.Q1y * T1y_1y() * pol2.U1y;
            U += pol1.Q1y * T1y_1z() * pol2.U1z;
            U += pol1.Q1z * T1z_1x() * pol2.U1x;
            U += pol1.Q1z * T1z_1y() * pol2.U1y;
            U += pol1.Q1z * T1z_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 0) {
        if (a*u3 < 40) {
            U += pol1.U1x * TU1x_1x() * pol2.Q1x;
            U += pol1.U1x * TU1x_1y() * pol2.Q1y;
            U += pol1.U1x * TU1x_1z() * pol2.Q1z;
            U += pol1.U1y * TU1y_1x() * pol2.Q1x;
            U += pol1.U1y * TU1y_1y() * pol2.Q1y;
            U += pol1.U1y * TU1y_1z() * pol2.Q1z;
            U += pol1.U1z * TU1z_1x() * pol2.Q1x;
            U += pol1.U1z * TU1z_1y() * pol2.Q1y;
            U += pol1.U1z * TU1z_1z() * pol2.Q1z;
        }
        else {
            U += pol1.U1x * T1x_1x() * pol2.Q1x;
            U += pol1.U1x * T1x_1y() * pol2.Q1y;
            U += pol1.U1x * T1x_1z() * pol2.Q1z;
            U += pol1.U1y * T1y_1x() * pol2.Q1x;
            U += pol1.U1y * T1y_1y() * pol2.Q1y;
            U += pol1.U1y * T1y_1z() * pol2.Q1z;
            U += pol1.U1z * T1z_1x() * pol2.Q1x;
            U += pol1.U1z * T1z_1y() * pol2.Q1y;
            U += pol1.U1z * T1z_1z() * pol2.Q1z;
        }
    }

    if (pol1._rank > 1 && pol2._rank > 0) {
        E += pol1.Q20 * T20_1x() * pol2.Q1x;
        E += pol1.Q20 * T20_1y() * pol2.Q1y;
        E += pol1.Q20 * T20_1z() * pol2.Q1z;

        E += pol1.Q21c * T21c_1x() * pol2.Q1x;
        E += pol1.Q21c * T21c_1y() * pol2.Q1y;
        E += pol1.Q21c * T21c_1z() * pol2.Q1z;

        E += pol1.Q21s * T21s_1x() * pol2.Q1x;
        E += pol1.Q21s * T21s_1y() * pol2.Q1y;
        E += pol1.Q21s * T21s_1z() * pol2.Q1z;

        E += pol1.Q22c * T22c_1x() * pol2.Q1x;
        E += pol1.Q22c * T22c_1y() * pol2.Q1y;
        E += pol1.Q22c * T22c_1z() * pol2.Q1z;

        E += pol1.Q22s * T22s_1x() * pol2.Q1x;
        E += pol1.Q22s * T22s_1y() * pol2.Q1y;
        E += pol1.Q22s * T22s_1z() * pol2.Q1z;
    }

    if (pol1._rank > 0 && pol2._rank > 1) {
        E += pol1.Q1x * T1x_20() * pol2.Q20;
        E += pol1.Q1y * T1y_20() * pol2.Q20;
        E += pol1.Q1z * T1z_20() * pol2.Q20;

        E += pol1.Q1x * T1x_21c() * pol2.Q21c;
        E += pol1.Q1y * T1y_21c() * pol2.Q21c;
        E += pol1.Q1z * T1z_21c() * pol2.Q21c;

        E += pol1.Q1x * T1x_21s() * pol2.Q21s;
        E += pol1.Q1y * T1y_21s() * pol2.Q21s;
        E += pol1.Q1z * T1z_21s() * pol2.Q21s;

        E += pol1.Q1x * T1x_22c() * pol2.Q22c;
        E += pol1.Q1y * T1y_22c() * pol2.Q22c;
        E += pol1.Q1z * T1z_22c() * pol2.Q22c;

        E += pol1.Q1x * T1x_22s() * pol2.Q22s;
        E += pol1.Q1y * T1y_22s() * pol2.Q22s;
        E += pol1.Q1z * T1z_22s() * pol2.Q22s;
    }

    if (pol1._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.Q20  * TU20_1x()  * pol2.U1x;
            U += pol1.Q20  * TU20_1y()  * pol2.U1y;
            U += pol1.Q20  * TU20_1z()  * pol2.U1z;
            U += pol1.Q21c * TU21c_1x() * pol2.U1x;
            U += pol1.Q21c * TU21c_1y() * pol2.U1y;
            U += pol1.Q21c * TU21c_1z() * pol2.U1z;
            U += pol1.Q21s * TU21s_1x() * pol2.U1x;
            U += pol1.Q21s * TU21s_1y() * pol2.U1y;
            U += pol1.Q21s * TU21s_1z() * pol2.U1z;
            U += pol1.Q22c * TU22c_1x() * pol2.U1x;
            U += pol1.Q22c * TU22c_1y() * pol2.U1y;
            U += pol1.Q22c * TU22c_1z() * pol2.U1z;
            U += pol1.Q22s * TU22s_1x() * pol2.U1x;
            U += pol1.Q22s * TU22s_1y() * pol2.U1y;
            U += pol1.Q22s * TU22s_1z() * pol2.U1z;
        }
        else {
            U += pol1.Q20  * T20_1x()  * pol2.U1x;
            U += pol1.Q20  * T20_1y()  * pol2.U1y;
            U += pol1.Q20  * T20_1z()  * pol2.U1z;
            U += pol1.Q21c * T21c_1x() * pol2.U1x;
            U += pol1.Q21c * T21c_1y() * pol2.U1y;
            U += pol1.Q21c * T21c_1z() * pol2.U1z;
            U += pol1.Q21s * T21s_1x() * pol2.U1x;
            U += pol1.Q21s * T21s_1y() * pol2.U1y;
            U += pol1.Q21s * T21s_1z() * pol2.U1z;
            U += pol1.Q22c * T22c_1x() * pol2.U1x;
            U += pol1.Q22c * T22c_1y() * pol2.U1y;
            U += pol1.Q22c * T22c_1z() * pol2.U1z;
            U += pol1.Q22s * T22s_1x() * pol2.U1x;
            U += pol1.Q22s * T22s_1y() * pol2.U1y;
            U += pol1.Q22s * T22s_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 1) {
        if (a*u3 < 40.0) {
            U += pol1.U1x * TU1x_20()  * pol2.Q20;
            U += pol1.U1x * TU1x_21c() * pol2.Q21c;
            U += pol1.U1x * TU1x_21s() * pol2.Q21s;
            U += pol1.U1x * TU1x_22c() * pol2.Q22c;
            U += pol1.U1x * TU1x_22s() * pol2.Q22s;
            U += pol1.U1y * TU1y_20()  * pol2.Q20;
            U += pol1.U1y * TU1y_21c() * pol2.Q21c;
            U += pol1.U1y * TU1y_21s() * pol2.Q21s;
            U += pol1.U1y * TU1y_22c() * pol2.Q22c;
            U += pol1.U1y * TU1y_22s() * pol2.Q22s;
            U += pol1.U1z * TU1z_20()  * pol2.Q20;
            U += pol1.U1z * TU1z_21c() * pol2.Q21c;
            U += pol1.U1z * TU1z_21s() * pol2.Q21s;
            U += pol1.U1z * TU1z_22c() * pol2.Q22c;
            U += pol1.U1z * TU1z_22s() * pol2.Q22s;
        }
        else {
            U += pol1.U1x * T1x_20()  * pol2.Q20;
            U += pol1.U1x * T1x_21c() * pol2.Q21c;
            U += pol1.U1x * T1x_21s() * pol2.Q21s;
            U += pol1.U1x * T1x_22c() * pol2.Q22c;
            U += pol1.U1x * T1x_22s() * pol2.Q22s;
            U += pol1.U1y * T1y_20()  * pol2.Q20;
            U += pol1.U1y * T1y_21c() * pol2.Q21c;
            U += pol1.U1y * T1y_21s() * pol2.Q21s;
            U += pol1.U1y * T1y_22c() * pol2.Q22c;
            U += pol1.U1y * T1y_22s() * pol2.Q22s;
            U += pol1.U1z * T1z_20()  * pol2.Q20;
            U += pol1.U1z * T1z_21c() * pol2.Q21c;
            U += pol1.U1z * T1z_21s() * pol2.Q21s;
            U += pol1.U1z * T1z_22c() * pol2.Q22c;
            U += pol1.U1z * T1z_22s() * pol2.Q22s;
        }
    }

    if (pol1._rank > 1 && pol2._rank > 1) {
        E += pol1.Q20  * T20_20()   * pol2.Q20;
        E += pol1.Q21c * T21c_21c() * pol2.Q21c;
        E += pol1.Q21s * T21s_21s() * pol2.Q21s;
        E += pol1.Q22c * T22c_22c() * pol2.Q22c;
        E += pol1.Q22s * T22s_22s() * pol2.Q22s;


        E += pol1.Q20  * T20_21c() * pol2.Q21c;
        E += pol1.Q20  * T20_21s() * pol2.Q21s;
        E += pol1.Q20  * T20_22c() * pol2.Q22c;
        E += pol1.Q20  * T20_22s() * pol2.Q22s;
        E += pol1.Q21c * T21c_20() * pol2.Q20;
        E += pol1.Q21s * T21s_20() * pol2.Q20;
        E += pol1.Q22c * T22c_20() * pol2.Q20;
        E += pol1.Q22s * T22s_20() * pol2.Q20;


        E += pol1.Q21c * T21c_21s() * pol2.Q21s;
        E += pol1.Q21c * T21c_22c() * pol2.Q22c;
        E += pol1.Q21c * T21c_22s() * pol2.Q22s;
        E += pol1.Q21s * T21s_21c() * pol2.Q21c;
        E += pol1.Q22c * T22c_21c() * pol2.Q21c;
        E += pol1.Q22s * T22s_21c() * pol2.Q21c;


        E += pol1.Q21s * T21s_22c() * pol2.Q22c;
        E += pol1.Q21s * T21s_22s() * pol2.Q22s;
        E += pol1.Q22c * T22c_21s() * pol2.Q21s;
        E += pol1.Q22s * T22s_21s() * pol2.Q21s;

        E += pol1.Q22s * T22s_22c() * pol2.Q22c;
        E += pol1.Q22c * T22c_22s() * pol2.Q22s;
    }

    // Take into account work required to induce multipoles
    U *= 0.5;

    return E + U;
}


inline double XInteractor::E_f(APolarSite &pol1, APolarSite &pol2) {

//    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
//    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
//    //e12  = pol2.getPos() - pol1.getPos();
//    R    = 1/abs(e12);
//    R2   = R*R;
//    R3   = R2*R;
//    R4   = R3*R;
//    R5   = R4*R;
//    e12 *= R;
//
//    // Thole damping init.
//    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));
//
//
//    //cout << "frag1 " << pol1.getFragment()->getId() << endl;
//    //cout << "frag2 " << pol2.getFragment()->getId() << endl;
//    //cout << "seg1  " << pol1.getSegment()->getId() << endl;
//    //cout << "seg2  " << pol2.getSegment()->getId() << endl;
//
//
////        rax =   pol1._locX * e12;
////        ray =   pol1._locY * e12;
////        raz =   pol1._locZ * e12;
////        rbx = - pol2._locX * e12;
////        rby = - pol2._locY * e12;
////        rbz = - pol2._locZ * e12;
//
//        rax = e12.getX();
//        ray = e12.getY();
//        raz = e12.getZ();
//        rbx = - rax;
//        rby = - ray;
//        rbz = - raz;
//
////        cxx = pol1._locX * pol2._locX;
////        cxy = pol1._locX * pol2._locY;
////        cxz = pol1._locX * pol2._locZ;
////        cyx = pol1._locY * pol2._locX;
////        cyy = pol1._locY * pol2._locY;
////        cyz = pol1._locY * pol2._locZ;
////        czx = pol1._locZ * pol2._locX;
////        czy = pol1._locZ * pol2._locY;
////        czz = pol1._locZ * pol2._locZ;
//
//        cxx = 1;
//        cxy = 0;
//        cxz = 0;
//        cyx = 0;
//        cyy = 1;
//        cyz = 0;
//        czx = 0;
//        czy = 0;
//        czz = 1;

    double epp = 0.0; // <- Interaction perm. <> perm.
    double epu = 0.0; // <- Interaction perm. <> induced
    double euu = 0.0; // <- Interaction induced <> induced>
    
    
        //cout << "r1  " << pol1.getPos() << endl;
        //cout << "r2  " << pol2.getPos() << endl;
        //cout << "R   " << 1/R << endl;
        //cout << "e12 " << e12 << endl;

        epp += pol1.Q00 * T00_00() * pol2.Q00;

        //cout << "E up to q <-> q " << E << endl;

    if (a*u3 < 40) {
        epu += pol1.U1x * TU1x_00() * pol2.Q00;
        epu += pol1.U1y * TU1y_00() * pol2.Q00;
        epu += pol1.U1z * TU1z_00() * pol2.Q00;

        epu += pol1.Q00 * TU00_1x() * pol2.U1x;
        epu += pol1.Q00 * TU00_1y() * pol2.U1y;
        epu += pol1.Q00 * TU00_1z() * pol2.U1z;
    }
    else {
        epu += pol1.U1x * T1x_00() * pol2.Q00;
        epu += pol1.U1y * T1y_00() * pol2.Q00;
        epu += pol1.U1z * T1z_00() * pol2.Q00;

        epu += pol1.Q00 * T00_1x() * pol2.U1x;
        epu += pol1.Q00 * T00_1y() * pol2.U1y;
        epu += pol1.Q00 * T00_1z() * pol2.U1z;
    }



    if (pol1._rank > 0) {
        epp += pol1.Q1x * T1x_00() * pol2.Q00;
        //cout << "E1x_00 " << pol1.Q1x * T1x_00() * pol2.Q00 << endl;
        epp += pol1.Q1y * T1y_00() * pol2.Q00;
        //cout << "E1y_00 " << pol1.Q1y * T1y_00() * pol2.Q00 << endl;
        epp += pol1.Q1z * T1z_00() * pol2.Q00;
        //cout << "E1z_00 " << pol1.Q1z * T1z_00() * pol2.Q00 << endl;
    }

    if (pol2._rank > 0) {
        epp += pol1.Q00 * T00_1x() * pol2.Q1x;
        //cout << "E00_1x " << pol1.Q00 * T00_1x() * pol2.Q1x << endl;
        epp += pol1.Q00 * T00_1y() * pol2.Q1y;
        //cout << "E00_1y " << pol1.Q00 * T00_1y() * pol2.Q1y << endl;
        epp += pol1.Q00 * T00_1z() * pol2.Q1z;
        //cout << "E00_1z " << pol1.Q00 * T00_1z() * pol2.Q1z << endl;
    }
        //cout << "E up to q <-> d " << E << endl;

    if (pol1._rank > 1) {
        epp += pol1.Q20  * T20_00()  * pol2.Q00;
        epp += pol1.Q21c * T21c_00() * pol2.Q00;
        epp += pol1.Q21s * T21s_00() * pol2.Q00;
        epp += pol1.Q22c * T22c_00() * pol2.Q00;
        epp += pol1.Q22s * T22s_00() * pol2.Q00;
    }

    if (pol2._rank > 1) {
        epp += pol1.Q00 * T00_20()  * pol2.Q20;
        epp += pol1.Q00 * T00_21c() * pol2.Q21c;
        epp += pol1.Q00 * T00_21s() * pol2.Q21s;
        epp += pol1.Q00 * T00_22c() * pol2.Q22c;
        epp += pol1.Q00 * T00_22s() * pol2.Q22s;
    }
        //cout << "E up to q <-> Q " << E << endl;

    if (pol1._rank > 0 && pol2._rank > 0) {
        epp += pol1.Q1x * T1x_1x() * pol2.Q1x;
        //cout << "E1x_1x " << pol1.Q1x * T1x_1x() * pol2.Q1x << endl;
        epp += pol1.Q1x * T1x_1y() * pol2.Q1y;
        //cout << "E1x_1y " << pol1.Q1x * T1x_1y() * pol2.Q1y << endl;
        epp += pol1.Q1x * T1x_1z() * pol2.Q1z;
        //cout << "E1x_1z " << pol1.Q1x * T1x_1z() * pol2.Q1z << endl;

        epp += pol1.Q1y * T1y_1x() * pol2.Q1x;
        //cout << "E1y_1x " << pol1.Q1y * T1y_1x() * pol2.Q1x << endl;
        epp += pol1.Q1y * T1y_1y() * pol2.Q1y;
        //cout << "E1y_1y " << pol1.Q1y * T1y_1y() * pol2.Q1y << endl;
        epp += pol1.Q1y * T1y_1z() * pol2.Q1z;
        //cout << "E1y_1z " << pol1.Q1y * T1y_1z() * pol2.Q1z << endl;

        epp += pol1.Q1z * T1z_1x() * pol2.Q1x;
        //cout << "E1z_1x " << pol1.Q1z * T1z_1x() * pol2.Q1x << endl;
        epp += pol1.Q1z * T1z_1y() * pol2.Q1y;
        //cout << "E1z_1y " << pol1.Q1z * T1z_1y() * pol2.Q1y << endl;
        epp += pol1.Q1z * T1z_1z() * pol2.Q1z;
        //cout << "E1z_1z " << pol1.Q1z * T1z_1z() * pol2.Q1z << endl;
    }

    if (pol1._rank > 0) {
        if (a*u3 < 40) {
            epu += pol1.Q1x * TU1x_1x() * pol2.U1x;
            epu += pol1.Q1x * TU1x_1y() * pol2.U1y;
            epu += pol1.Q1x * TU1x_1z() * pol2.U1z;
            epu += pol1.Q1y * TU1y_1x() * pol2.U1x;
            epu += pol1.Q1y * TU1y_1y() * pol2.U1y;
            epu += pol1.Q1y * TU1y_1z() * pol2.U1z;
            epu += pol1.Q1z * TU1z_1x() * pol2.U1x;
            epu += pol1.Q1z * TU1z_1y() * pol2.U1y;
            epu += pol1.Q1z * TU1z_1z() * pol2.U1z;
        }
        else {
            epu += pol1.Q1x * T1x_1x() * pol2.U1x;
            epu += pol1.Q1x * T1x_1y() * pol2.U1y;
            epu += pol1.Q1x * T1x_1z() * pol2.U1z;
            epu += pol1.Q1y * T1y_1x() * pol2.U1x;
            epu += pol1.Q1y * T1y_1y() * pol2.U1y;
            epu += pol1.Q1y * T1y_1z() * pol2.U1z;
            epu += pol1.Q1z * T1z_1x() * pol2.U1x;
            epu += pol1.Q1z * T1z_1y() * pol2.U1y;
            epu += pol1.Q1z * T1z_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 0) {
        if (a*u3 < 40) {
            epu += pol1.U1x * TU1x_1x() * pol2.Q1x;
            epu += pol1.U1x * TU1x_1y() * pol2.Q1y;
            epu += pol1.U1x * TU1x_1z() * pol2.Q1z;
            epu += pol1.U1y * TU1y_1x() * pol2.Q1x;
            epu += pol1.U1y * TU1y_1y() * pol2.Q1y;
            epu += pol1.U1y * TU1y_1z() * pol2.Q1z;
            epu += pol1.U1z * TU1z_1x() * pol2.Q1x;
            epu += pol1.U1z * TU1z_1y() * pol2.Q1y;
            epu += pol1.U1z * TU1z_1z() * pol2.Q1z;
        }
        else {
            epu += pol1.U1x * T1x_1x() * pol2.Q1x;
            epu += pol1.U1x * T1x_1y() * pol2.Q1y;
            epu += pol1.U1x * T1x_1z() * pol2.Q1z;
            epu += pol1.U1y * T1y_1x() * pol2.Q1x;
            epu += pol1.U1y * T1y_1y() * pol2.Q1y;
            epu += pol1.U1y * T1y_1z() * pol2.Q1z;
            epu += pol1.U1z * T1z_1x() * pol2.Q1x;
            epu += pol1.U1z * T1z_1y() * pol2.Q1y;
            epu += pol1.U1z * T1z_1z() * pol2.Q1z;
        }
    }
        //cout << "E up to d <-> d " << E << endl;

    if (pol1._rank > 1 && pol2._rank > 0) {
        epp += pol1.Q20 * T20_1x() * pol2.Q1x;
        epp += pol1.Q20 * T20_1y() * pol2.Q1y;
        epp += pol1.Q20 * T20_1z() * pol2.Q1z;

        epp += pol1.Q21c * T21c_1x() * pol2.Q1x;
        epp += pol1.Q21c * T21c_1y() * pol2.Q1y;
        epp += pol1.Q21c * T21c_1z() * pol2.Q1z;

        epp += pol1.Q21s * T21s_1x() * pol2.Q1x;
        epp += pol1.Q21s * T21s_1y() * pol2.Q1y;
        epp += pol1.Q21s * T21s_1z() * pol2.Q1z;

        epp += pol1.Q22c * T22c_1x() * pol2.Q1x;
        epp += pol1.Q22c * T22c_1y() * pol2.Q1y;
        epp += pol1.Q22c * T22c_1z() * pol2.Q1z;

        epp += pol1.Q22s * T22s_1x() * pol2.Q1x;
        epp += pol1.Q22s * T22s_1y() * pol2.Q1y;
        epp += pol1.Q22s * T22s_1z() * pol2.Q1z;
    }

    if (pol1._rank > 0 && pol2._rank > 1) {
        epp += pol1.Q1x * T1x_20() * pol2.Q20;
        epp += pol1.Q1y * T1y_20() * pol2.Q20;
        epp += pol1.Q1z * T1z_20() * pol2.Q20;

        epp += pol1.Q1x * T1x_21c() * pol2.Q21c;
        epp += pol1.Q1y * T1y_21c() * pol2.Q21c;
        epp += pol1.Q1z * T1z_21c() * pol2.Q21c;

        epp += pol1.Q1x * T1x_21s() * pol2.Q21s;
        epp += pol1.Q1y * T1y_21s() * pol2.Q21s;
        epp += pol1.Q1z * T1z_21s() * pol2.Q21s;

        epp += pol1.Q1x * T1x_22c() * pol2.Q22c;
        epp += pol1.Q1y * T1y_22c() * pol2.Q22c;
        epp += pol1.Q1z * T1z_22c() * pol2.Q22c;

        epp += pol1.Q1x * T1x_22s() * pol2.Q22s;
        epp += pol1.Q1y * T1y_22s() * pol2.Q22s;
        epp += pol1.Q1z * T1z_22s() * pol2.Q22s;
    }

    if (pol1._rank > 1) {
        if (a*u3 < 40.0) {
            epu += pol1.Q20  * TU20_1x()  * pol2.U1x;
            epu += pol1.Q20  * TU20_1y()  * pol2.U1y;
            epu += pol1.Q20  * TU20_1z()  * pol2.U1z;
            epu += pol1.Q21c * TU21c_1x() * pol2.U1x;
            epu += pol1.Q21c * TU21c_1y() * pol2.U1y;
            epu += pol1.Q21c * TU21c_1z() * pol2.U1z;
            epu += pol1.Q21s * TU21s_1x() * pol2.U1x;
            epu += pol1.Q21s * TU21s_1y() * pol2.U1y;
            epu += pol1.Q21s * TU21s_1z() * pol2.U1z;
            epu += pol1.Q22c * TU22c_1x() * pol2.U1x;
            epu += pol1.Q22c * TU22c_1y() * pol2.U1y;
            epu += pol1.Q22c * TU22c_1z() * pol2.U1z;
            epu += pol1.Q22s * TU22s_1x() * pol2.U1x;
            epu += pol1.Q22s * TU22s_1y() * pol2.U1y;
            epu += pol1.Q22s * TU22s_1z() * pol2.U1z;
        }
        else {
            epu += pol1.Q20  * T20_1x()  * pol2.U1x;
            epu += pol1.Q20  * T20_1y()  * pol2.U1y;
            epu += pol1.Q20  * T20_1z()  * pol2.U1z;
            epu += pol1.Q21c * T21c_1x() * pol2.U1x;
            epu += pol1.Q21c * T21c_1y() * pol2.U1y;
            epu += pol1.Q21c * T21c_1z() * pol2.U1z;
            epu += pol1.Q21s * T21s_1x() * pol2.U1x;
            epu += pol1.Q21s * T21s_1y() * pol2.U1y;
            epu += pol1.Q21s * T21s_1z() * pol2.U1z;
            epu += pol1.Q22c * T22c_1x() * pol2.U1x;
            epu += pol1.Q22c * T22c_1y() * pol2.U1y;
            epu += pol1.Q22c * T22c_1z() * pol2.U1z;
            epu += pol1.Q22s * T22s_1x() * pol2.U1x;
            epu += pol1.Q22s * T22s_1y() * pol2.U1y;
            epu += pol1.Q22s * T22s_1z() * pol2.U1z;
        }
    }
    if (pol2._rank > 1) {
        if (a*u3 < 40.0) {
            epu += pol1.U1x * TU1x_20()  * pol2.Q20;
            epu += pol1.U1x * TU1x_21c() * pol2.Q21c;
            epu += pol1.U1x * TU1x_21s() * pol2.Q21s;
            epu += pol1.U1x * TU1x_22c() * pol2.Q22c;
            epu += pol1.U1x * TU1x_22s() * pol2.Q22s;
            epu += pol1.U1y * TU1y_20()  * pol2.Q20;
            epu += pol1.U1y * TU1y_21c() * pol2.Q21c;
            epu += pol1.U1y * TU1y_21s() * pol2.Q21s;
            epu += pol1.U1y * TU1y_22c() * pol2.Q22c;
            epu += pol1.U1y * TU1y_22s() * pol2.Q22s;
            epu += pol1.U1z * TU1z_20()  * pol2.Q20;
            epu += pol1.U1z * TU1z_21c() * pol2.Q21c;
            epu += pol1.U1z * TU1z_21s() * pol2.Q21s;
            epu += pol1.U1z * TU1z_22c() * pol2.Q22c;
            epu += pol1.U1z * TU1z_22s() * pol2.Q22s;
        }
        else {
            epu += pol1.U1x * T1x_20()  * pol2.Q20;
            epu += pol1.U1x * T1x_21c() * pol2.Q21c;
            epu += pol1.U1x * T1x_21s() * pol2.Q21s;
            epu += pol1.U1x * T1x_22c() * pol2.Q22c;
            epu += pol1.U1x * T1x_22s() * pol2.Q22s;
            epu += pol1.U1y * T1y_20()  * pol2.Q20;
            epu += pol1.U1y * T1y_21c() * pol2.Q21c;
            epu += pol1.U1y * T1y_21s() * pol2.Q21s;
            epu += pol1.U1y * T1y_22c() * pol2.Q22c;
            epu += pol1.U1y * T1y_22s() * pol2.Q22s;
            epu += pol1.U1z * T1z_20()  * pol2.Q20;
            epu += pol1.U1z * T1z_21c() * pol2.Q21c;
            epu += pol1.U1z * T1z_21s() * pol2.Q21s;
            epu += pol1.U1z * T1z_22c() * pol2.Q22c;
            epu += pol1.U1z * T1z_22s() * pol2.Q22s;
        }
    }
        //cout << "E up to d <-> Q " << E << endl;

    if (pol1._rank > 1 && pol2._rank > 1) {
        epp += pol1.Q20  * T20_20()   * pol2.Q20;
        epp += pol1.Q21c * T21c_21c() * pol2.Q21c;
        epp += pol1.Q21s * T21s_21s() * pol2.Q21s;
        epp += pol1.Q22c * T22c_22c() * pol2.Q22c;
        epp += pol1.Q22s * T22s_22s() * pol2.Q22s;


        epp += pol1.Q20  * T20_21c() * pol2.Q21c;
        epp += pol1.Q20  * T20_21s() * pol2.Q21s;
        epp += pol1.Q20  * T20_22c() * pol2.Q22c;
        epp += pol1.Q20  * T20_22s() * pol2.Q22s;
        epp += pol1.Q21c * T21c_20() * pol2.Q20;
        epp += pol1.Q21s * T21s_20() * pol2.Q20;
        epp += pol1.Q22c * T22c_20() * pol2.Q20;
        epp += pol1.Q22s * T22s_20() * pol2.Q20;


        epp += pol1.Q21c * T21c_21s() * pol2.Q21s;
        epp += pol1.Q21c * T21c_22c() * pol2.Q22c;
        epp += pol1.Q21c * T21c_22s() * pol2.Q22s;
        epp += pol1.Q21s * T21s_21c() * pol2.Q21c;
        epp += pol1.Q22c * T22c_21c() * pol2.Q21c;
        epp += pol1.Q22s * T22s_21c() * pol2.Q21c;


        epp += pol1.Q21s * T21s_22c() * pol2.Q22c;
        epp += pol1.Q21s * T21s_22s() * pol2.Q22s;
        epp += pol1.Q22c * T22c_21s() * pol2.Q21s;
        epp += pol1.Q22s * T22s_21s() * pol2.Q21s;

        epp += pol1.Q22s * T22s_22c() * pol2.Q22c;
        epp += pol1.Q22c * T22c_22s() * pol2.Q22s;
    }
        //cout << "E up to Q <-> Q " << E << endl;      
        
    
        
    // Induced <> induced interaction
    // ... Note that this contribution is canceled by the induction work.        
    if (a*u3 < 40) {
        euu += pol1.U1x * TU1x_1x() * pol2.U1x;
        euu += pol1.U1x * TU1x_1y() * pol2.U1y;
        euu += pol1.U1x * TU1x_1z() * pol2.U1z;
        euu += pol1.U1y * TU1y_1x() * pol2.U1x;
        euu += pol1.U1y * TU1y_1y() * pol2.U1y;
        euu += pol1.U1y * TU1y_1z() * pol2.U1z;
        euu += pol1.U1z * TU1z_1x() * pol2.U1x;
        euu += pol1.U1z * TU1z_1y() * pol2.U1y;
        euu += pol1.U1z * TU1z_1z() * pol2.U1z;
    }
    else {
        euu += pol1.U1x * T1x_1x() * pol2.U1x;
        euu += pol1.U1x * T1x_1y() * pol2.U1y;
        euu += pol1.U1x * T1x_1z() * pol2.U1z;
        euu += pol1.U1y * T1y_1x() * pol2.U1x;
        euu += pol1.U1y * T1y_1y() * pol2.U1y;
        euu += pol1.U1y * T1y_1z() * pol2.U1z;
        euu += pol1.U1z * T1z_1x() * pol2.U1x;
        euu += pol1.U1z * T1z_1y() * pol2.U1y;
        euu += pol1.U1z * T1z_1z() * pol2.U1z;
    }
    
    EP += epp;
    EU_INTER += 0.5 * epu;
    
    EPP += epp;
    EPU += epu;
    EUU += euu;
    
    return epp + epu + euu;
}


inline double XInteractor::E_m(APolarSite &pol1, APolarSite &pol2) {
    
//    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
//    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
//    //e12  = pol2.getPos() - pol1.getPos();
//    R    = 1/abs(e12);
//    R2   = R*R;
//    R3   = R2*R;
//    R4   = R3*R;
//    R5   = R4*R;
//    e12 *= R;
//
//    // Thole damping init.
//    u3   = 1 / (R3 * sqrt(pol1.eigendamp * pol2.eigendamp));
//
//
//    //cout << "frag1 " << pol1.getFragment()->getId() << endl;
//    //cout << "frag2 " << pol2.getFragment()->getId() << endl;
//    //cout << "seg1  " << pol1.getSegment()->getId() << endl;
//    //cout << "seg2  " << pol2.getSegment()->getId() << endl;
//
//
////        rax =   pol1._locX * e12;
////        ray =   pol1._locY * e12;
////        raz =   pol1._locZ * e12;
////        rbx = - pol2._locX * e12;
////        rby = - pol2._locY * e12;
////        rbz = - pol2._locZ * e12;
//
//        rax = e12.getX();
//        ray = e12.getY();
//        raz = e12.getZ();
//        rbx = - rax;
//        rby = - ray;
//        rbz = - raz;
//
////        cxx = pol1._locX * pol2._locX;
////        cxy = pol1._locX * pol2._locY;
////        cxz = pol1._locX * pol2._locZ;
////        cyx = pol1._locY * pol2._locX;
////        cyy = pol1._locY * pol2._locY;
////        cyz = pol1._locY * pol2._locZ;
////        czx = pol1._locZ * pol2._locX;
////        czy = pol1._locZ * pol2._locY;
////        czz = pol1._locZ * pol2._locZ;
//
//        cxx = 1;
//        cxy = 0;
//        cxz = 0;
//        cyx = 0;
//        cyy = 1;
//        cyz = 0;
//        czx = 0;
//        czy = 0;
//        czz = 1;

        
    double epu = 0.0; // <- Interaction perm. <> induced
    double euu = 0.0; // <- Interaction induced <> induced>


    if (a*u3 < 40) {
        epu += pol1.U1x * TU1x_00() * pol2.Q00;
        epu += pol1.U1y * TU1y_00() * pol2.Q00;
        epu += pol1.U1z * TU1z_00() * pol2.Q00;

    }
    else {
        epu += pol1.U1x * T1x_00() * pol2.Q00;
        epu += pol1.U1y * T1y_00() * pol2.Q00;
        epu += pol1.U1z * T1z_00() * pol2.Q00;
    }

    if (pol2._rank > 0) {
        if (a*u3 < 40) {
            epu += pol1.U1x * TU1x_1x() * pol2.Q1x;
            epu += pol1.U1x * TU1x_1y() * pol2.Q1y;
            epu += pol1.U1x * TU1x_1z() * pol2.Q1z;
            epu += pol1.U1y * TU1y_1x() * pol2.Q1x;
            epu += pol1.U1y * TU1y_1y() * pol2.Q1y;
            epu += pol1.U1y * TU1y_1z() * pol2.Q1z;
            epu += pol1.U1z * TU1z_1x() * pol2.Q1x;
            epu += pol1.U1z * TU1z_1y() * pol2.Q1y;
            epu += pol1.U1z * TU1z_1z() * pol2.Q1z;
        }
        else {
            epu += pol1.U1x * T1x_1x() * pol2.Q1x;
            epu += pol1.U1x * T1x_1y() * pol2.Q1y;
            epu += pol1.U1x * T1x_1z() * pol2.Q1z;
            epu += pol1.U1y * T1y_1x() * pol2.Q1x;
            epu += pol1.U1y * T1y_1y() * pol2.Q1y;
            epu += pol1.U1y * T1y_1z() * pol2.Q1z;
            epu += pol1.U1z * T1z_1x() * pol2.Q1x;
            epu += pol1.U1z * T1z_1y() * pol2.Q1y;
            epu += pol1.U1z * T1z_1z() * pol2.Q1z;
        }
    }

    if (pol2._rank > 1) {
        if (a*u3 < 40.0) {
            epu += pol1.U1x * TU1x_20()  * pol2.Q20;
            epu += pol1.U1x * TU1x_21c() * pol2.Q21c;
            epu += pol1.U1x * TU1x_21s() * pol2.Q21s;
            epu += pol1.U1x * TU1x_22c() * pol2.Q22c;
            epu += pol1.U1x * TU1x_22s() * pol2.Q22s;
            epu += pol1.U1y * TU1y_20()  * pol2.Q20;
            epu += pol1.U1y * TU1y_21c() * pol2.Q21c;
            epu += pol1.U1y * TU1y_21s() * pol2.Q21s;
            epu += pol1.U1y * TU1y_22c() * pol2.Q22c;
            epu += pol1.U1y * TU1y_22s() * pol2.Q22s;
            epu += pol1.U1z * TU1z_20()  * pol2.Q20;
            epu += pol1.U1z * TU1z_21c() * pol2.Q21c;
            epu += pol1.U1z * TU1z_21s() * pol2.Q21s;
            epu += pol1.U1z * TU1z_22c() * pol2.Q22c;
            epu += pol1.U1z * TU1z_22s() * pol2.Q22s;
        }
        else {
            epu += pol1.U1x * T1x_20()  * pol2.Q20;
            epu += pol1.U1x * T1x_21c() * pol2.Q21c;
            epu += pol1.U1x * T1x_21s() * pol2.Q21s;
            epu += pol1.U1x * T1x_22c() * pol2.Q22c;
            epu += pol1.U1x * T1x_22s() * pol2.Q22s;
            epu += pol1.U1y * T1y_20()  * pol2.Q20;
            epu += pol1.U1y * T1y_21c() * pol2.Q21c;
            epu += pol1.U1y * T1y_21s() * pol2.Q21s;
            epu += pol1.U1y * T1y_22c() * pol2.Q22c;
            epu += pol1.U1y * T1y_22s() * pol2.Q22s;
            epu += pol1.U1z * T1z_20()  * pol2.Q20;
            epu += pol1.U1z * T1z_21c() * pol2.Q21c;
            epu += pol1.U1z * T1z_21s() * pol2.Q21s;
            epu += pol1.U1z * T1z_22c() * pol2.Q22c;
            epu += pol1.U1z * T1z_22s() * pol2.Q22s;
        }
    }
               
    if (a*u3 < 40) {
        euu += pol1.U1x * TU1x_1x() * pol2.U1x;
        euu += pol1.U1x * TU1x_1y() * pol2.U1y;
        euu += pol1.U1x * TU1x_1z() * pol2.U1z;
        euu += pol1.U1y * TU1y_1x() * pol2.U1x;
        euu += pol1.U1y * TU1y_1y() * pol2.U1y;
        euu += pol1.U1y * TU1y_1z() * pol2.U1z;
        euu += pol1.U1z * TU1z_1x() * pol2.U1x;
        euu += pol1.U1z * TU1z_1y() * pol2.U1y;
        euu += pol1.U1z * TU1z_1z() * pol2.U1z;
    }
    else {
        euu += pol1.U1x * T1x_1x() * pol2.U1x;
        euu += pol1.U1x * T1x_1y() * pol2.U1y;
        euu += pol1.U1x * T1x_1z() * pol2.U1z;
        euu += pol1.U1y * T1y_1x() * pol2.U1x;
        euu += pol1.U1y * T1y_1y() * pol2.U1y;
        euu += pol1.U1y * T1y_1z() * pol2.U1z;
        euu += pol1.U1z * T1z_1x() * pol2.U1x;
        euu += pol1.U1z * T1z_1y() * pol2.U1y;
        euu += pol1.U1z * T1z_1z() * pol2.U1z;
    }
    
    euu *= -0.5;
    epu *= -0.5;

    EPU += epu;
    EUU += euu;
    
    return epu + euu;
}


inline double XInteractor::E_QQ_ERFC(APolarSite &pol1, APolarSite &pol2, double &ew_alpha) {
    e12  = pol2.getPos() - pol1.getPos();    
    R    = 1/abs(e12);
    return pol1.Q00 * pol2.Q00 * R * erfc(ew_alpha/R);
    //return pol1.Q00 * T00_00() * pol2.Q00;
}


inline double XInteractor::E_QQ_ERF(APolarSite &pol1, APolarSite &pol2, double &ew_alpha) {
    e12  = pol2.getPos() - pol1.getPos();
    R = (abs(e12) < 1e-50) ? 1e+50 : 1/abs(e12);    
    return pol1.Q00 * pol2.Q00 * R * erf(ew_alpha/R);
}


inline double XInteractor::E_QQ_K0(APolarSite &pol1, APolarSite &pol2, double &ew_alpha) {
    // NOTE Without prefactor 2*PI/(Lx*Ly)
    e12  = pol2.getPos() - pol1.getPos();
    double z12  = e12.getZ();
    return - pol1.Q00 * pol2.Q00 * (   exp(-ew_alpha*ew_alpha*z12*z12) / (sqrt(M_PI)*ew_alpha)   +    z12*erf(ew_alpha*z12)    );  
}


inline double XInteractor::E_QQ_KK(APolarSite &pol1, APolarSite &pol2, double &ew_alpha, vec &k) {
    // NOTE Without prefactor 2*PI/(Lx*Ly)
    e12  = pol2.getPos() - pol1.getPos();
    double z12  = e12.getZ();
    double K = abs(k);
    return pol1.Q00 * pol2.Q00 * cos(k*e12) / K * (    exp(K*z12)*erfc(K/(2*ew_alpha)+ew_alpha*z12)   +    exp(-K*z12)*erfc(K/(2*ew_alpha)-ew_alpha*z12)   ); 
}


inline double XInteractor::E_Q0_DQ(APolarSite &pol1, APolarSite &pol2) {
    
    // Counts these interactions:
    // pol1(q)       <> pol2(d,Q,...)
    // pol1(d,Q,...) <> pol2(q)
    
    double epp = 0.0;
    
    // d <> q
    if (pol1._rank > 0) {
        epp += pol1.Q1x * T1x_00() * pol2.Q00;
        epp += pol1.Q1y * T1y_00() * pol2.Q00;
        epp += pol1.Q1z * T1z_00() * pol2.Q00;
    }

    // q <> d
    if (pol2._rank > 0) {
        epp += pol1.Q00 * T00_1x() * pol2.Q1x;
        epp += pol1.Q00 * T00_1y() * pol2.Q1y;
        epp += pol1.Q00 * T00_1z() * pol2.Q1z;
    }

    // Q <> q
    if (pol1._rank > 1) {
        epp += pol1.Q20  * T20_00()  * pol2.Q00;
        epp += pol1.Q21c * T21c_00() * pol2.Q00;
        epp += pol1.Q21s * T21s_00() * pol2.Q00;
        epp += pol1.Q22c * T22c_00() * pol2.Q00;
        epp += pol1.Q22s * T22s_00() * pol2.Q00;
    }

    // q <> Q
    if (pol2._rank > 1) {
        epp += pol1.Q00 * T00_20()  * pol2.Q20;
        epp += pol1.Q00 * T00_21c() * pol2.Q21c;
        epp += pol1.Q00 * T00_21s() * pol2.Q21s;
        epp += pol1.Q00 * T00_22c() * pol2.Q22c;
        epp += pol1.Q00 * T00_22s() * pol2.Q22s;
    }

    // d <> d
    if (pol1._rank > 0 && pol2._rank > 0) {
        epp += pol1.Q1x * T1x_1x() * pol2.Q1x;
        epp += pol1.Q1x * T1x_1y() * pol2.Q1y;
        epp += pol1.Q1x * T1x_1z() * pol2.Q1z;

        epp += pol1.Q1y * T1y_1x() * pol2.Q1x;
        epp += pol1.Q1y * T1y_1y() * pol2.Q1y;
        epp += pol1.Q1y * T1y_1z() * pol2.Q1z;

        epp += pol1.Q1z * T1z_1x() * pol2.Q1x;
        epp += pol1.Q1z * T1z_1y() * pol2.Q1y;
        epp += pol1.Q1z * T1z_1z() * pol2.Q1z;
    }

    // Q <> d
    if (pol1._rank > 1 && pol2._rank > 0) {
        epp += pol1.Q20 * T20_1x() * pol2.Q1x;
        epp += pol1.Q20 * T20_1y() * pol2.Q1y;
        epp += pol1.Q20 * T20_1z() * pol2.Q1z;

        epp += pol1.Q21c * T21c_1x() * pol2.Q1x;
        epp += pol1.Q21c * T21c_1y() * pol2.Q1y;
        epp += pol1.Q21c * T21c_1z() * pol2.Q1z;

        epp += pol1.Q21s * T21s_1x() * pol2.Q1x;
        epp += pol1.Q21s * T21s_1y() * pol2.Q1y;
        epp += pol1.Q21s * T21s_1z() * pol2.Q1z;

        epp += pol1.Q22c * T22c_1x() * pol2.Q1x;
        epp += pol1.Q22c * T22c_1y() * pol2.Q1y;
        epp += pol1.Q22c * T22c_1z() * pol2.Q1z;

        epp += pol1.Q22s * T22s_1x() * pol2.Q1x;
        epp += pol1.Q22s * T22s_1y() * pol2.Q1y;
        epp += pol1.Q22s * T22s_1z() * pol2.Q1z;
    }

    // d <> Q
    if (pol1._rank > 0 && pol2._rank > 1) {
        epp += pol1.Q1x * T1x_20() * pol2.Q20;
        epp += pol1.Q1y * T1y_20() * pol2.Q20;
        epp += pol1.Q1z * T1z_20() * pol2.Q20;

        epp += pol1.Q1x * T1x_21c() * pol2.Q21c;
        epp += pol1.Q1y * T1y_21c() * pol2.Q21c;
        epp += pol1.Q1z * T1z_21c() * pol2.Q21c;

        epp += pol1.Q1x * T1x_21s() * pol2.Q21s;
        epp += pol1.Q1y * T1y_21s() * pol2.Q21s;
        epp += pol1.Q1z * T1z_21s() * pol2.Q21s;

        epp += pol1.Q1x * T1x_22c() * pol2.Q22c;
        epp += pol1.Q1y * T1y_22c() * pol2.Q22c;
        epp += pol1.Q1z * T1z_22c() * pol2.Q22c;

        epp += pol1.Q1x * T1x_22s() * pol2.Q22s;
        epp += pol1.Q1y * T1y_22s() * pol2.Q22s;
        epp += pol1.Q1z * T1z_22s() * pol2.Q22s;
    }

    // Q <> Q
    if (pol1._rank > 1 && pol2._rank > 1) {
        epp += pol1.Q20  * T20_20()   * pol2.Q20;
        epp += pol1.Q21c * T21c_21c() * pol2.Q21c;
        epp += pol1.Q21s * T21s_21s() * pol2.Q21s;
        epp += pol1.Q22c * T22c_22c() * pol2.Q22c;
        epp += pol1.Q22s * T22s_22s() * pol2.Q22s;


        epp += pol1.Q20  * T20_21c() * pol2.Q21c;
        epp += pol1.Q20  * T20_21s() * pol2.Q21s;
        epp += pol1.Q20  * T20_22c() * pol2.Q22c;
        epp += pol1.Q20  * T20_22s() * pol2.Q22s;
        epp += pol1.Q21c * T21c_20() * pol2.Q20;
        epp += pol1.Q21s * T21s_20() * pol2.Q20;
        epp += pol1.Q22c * T22c_20() * pol2.Q20;
        epp += pol1.Q22s * T22s_20() * pol2.Q20;


        epp += pol1.Q21c * T21c_21s() * pol2.Q21s;
        epp += pol1.Q21c * T21c_22c() * pol2.Q22c;
        epp += pol1.Q21c * T21c_22s() * pol2.Q22s;
        epp += pol1.Q21s * T21s_21c() * pol2.Q21c;
        epp += pol1.Q22c * T22c_21c() * pol2.Q21c;
        epp += pol1.Q22s * T22s_21c() * pol2.Q21c;


        epp += pol1.Q21s * T21s_22c() * pol2.Q22c;
        epp += pol1.Q21s * T21s_22s() * pol2.Q22s;
        epp += pol1.Q22c * T22c_21s() * pol2.Q21s;
        epp += pol1.Q22s * T22s_21s() * pol2.Q21s;

        epp += pol1.Q22s * T22s_22c() * pol2.Q22c;
        epp += pol1.Q22c * T22c_22s() * pol2.Q22s;
    }
    
    return epp;    
}


inline void XInteractor::BiasStat(APolarSite &pol1, APolarSite &pol2) {
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //          This implies that induced = - alpha * field
    //e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    e12  = pol2.getPos() - pol1.getPos();    
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

//        rax =   pol1._locX * e12;
//        ray =   pol1._locY * e12;
//        raz =   pol1._locZ * e12;
//        rbx = - pol2._locX * e12;
//        rby = - pol2._locY * e12;
//        rbz = - pol2._locZ * e12;

        rax = e12.getX(); rbx = - rax;
        ray = e12.getY(); rby = - ray;
        raz = e12.getZ(); rbz = - raz;

    if (pol1._rank > 0 || pol2._rank > 0) {
//        cxx = pol1._locX * pol2._locX;
//        cxy = pol1._locX * pol2._locY;
//        cxz = pol1._locX * pol2._locZ;
//        cyx = pol1._locY * pol2._locX;
//        cyy = pol1._locY * pol2._locY;
//        cyz = pol1._locY * pol2._locZ;
//        czx = pol1._locZ * pol2._locX;
//        czy = pol1._locZ * pol2._locY;
//        czz = pol1._locZ * pol2._locZ;

        cxx = 1;        cxy = 0;        cxz = 0;
        cyx = 0;        cyy = 1;        cyz = 0;
        czx = 0;        czy = 0;        czz = 1;
    }
}


inline void XInteractor::BiasIndu(APolarSite &pol1, APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    u3   = 1 / (R3 * sqrt( 
        1./3.*(pol1.Pxx*pol2.Pxx + pol1.Pxy*pol2.Pxy + pol1.Pxz*pol2.Pxz
             + pol1.Pxy*pol2.Pxy + pol1.Pyy*pol2.Pyy + pol1.Pyz*pol2.Pyz
             + pol1.Pxz*pol2.Pxz + pol1.Pyz*pol2.Pyz + pol1.Pzz*pol2.Pzz) ));

//        rax =   pol1._locX * e12;
//        ray =   pol1._locY * e12;
//        raz =   pol1._locZ * e12;
//        rbx = - pol2._locX * e12;
//        rby = - pol2._locY * e12;
//        rbz = - pol2._locZ * e12;

        rax = e12.getX(); rbx = - rax;
        ray = e12.getY(); rby = - ray;
        raz = e12.getZ(); rbz = - raz;

//        cxx = pol1._locX * pol2._locX;
//        cxy = pol1._locX * pol2._locY;
//        cxz = pol1._locX * pol2._locZ;
//        cyx = pol1._locY * pol2._locX;
//        cyy = pol1._locY * pol2._locY;
//        cyz = pol1._locY * pol2._locZ;
//        czx = pol1._locZ * pol2._locX;
//        czy = pol1._locZ * pol2._locY;
//        czz = pol1._locZ * pol2._locZ;

        cxx = 1;        cxy = 0;        cxz = 0;
        cyx = 0;        cyy = 1;        cyz = 0;
        czx = 0;        czy = 0;        czz = 1;

}


inline void XInteractor::RevBias() {
    
    e12 *= -1;
    
    rax *= -1;
    ray *= -1;
    raz *= -1;
    rbx *= -1;
    rby *= -1;
    rbz *= -1;
    
    // If C-matrix non-diagonal => transpose
    
}

}}

#endif

