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


#ifndef QMULTIPOLE_H
#define QMULTIPOLE_H


#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {


class QMultipole : public QMCalculator
{

public:

    QMultipole() {};
   ~QMultipole() {};

    string   Identify() { return "ZMultipole - Anisotropic Polarizabilities"; }

    // ++++++++++++++++++++++ //
    // Multipole Distribution //
    // ++++++++++++++++++++++ //

    void     Initialize(Property *options);
    void     EStatify(Topology *top, Property *options);
    void     DistributeMpoles(Topology *top);

    vector<APolarSite*> ParseGdmaFile(string filename, int state);
    vector< vec >      ParseCubeFileHeader(string filename);


    // ++++++++++++++++++++++++++++++++ //
    // Sub-Calculators ESP, ESF, MolPol //
    // ++++++++++++++++++++++++++++++++ //

    void     CalculateESP(vector< APolarSite* > &, FILE *);
    void     CalculateESPInput(Topology *top);
    void     CalculateESPRigid(Topology *top);
    void     CalculateESF(Topology *top);

    void     CalculateAlpha(Topology *top);
    void     CalculateAlphaInput(Topology *top);
    void     CalculateAlphaRigid(vector< APolarSite* > &, FILE *, bool);
    void     SiteAlphaInduce(Topology *, vector< APolarSite* > &);


    // ++++++++++++++++++++++++++++++++ //
    // Main Calculator => Site Operator //
    // ++++++++++++++++++++++++++++++++ //

    bool     EvaluateFrame(Topology *top);
    void     EvalSite(Topology *top, Segment *seg, int slot) { ; }
    void     Induce() { return; }
    void     CalcIntEnergy() { return; }

    void     InitSlotData(Topology *top) { ; }
    void     PostProcess(Topology *top) { ; }
    Segment *RequestNextSite(int opId, Topology *top);
    void     LockCout() { _coutMutex.Lock(); }
    void     UnlockCout() { _coutMutex.Unlock(); }


    // +++++++++++++++++++++++++++ //
    // Multipole Interaction Class //
    // +++++++++++++++++++++++++++ //

    class Interactor
    {
    public:

        Interactor(Topology *top, QMultipole *em) : _top(top), _em(em) {};
        Interactor() : _top(NULL), _em(NULL) {};
       ~Interactor() {};

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


        inline double EnergyInter(APolarSite &pol1, APolarSite &pol2);
        inline double EnergyInterESP(APolarSite &pol1, APolarSite &pol2);
        inline double EnergyIntra(APolarSite &pol1, APolarSite &pol2);
        inline void FieldPerm(APolarSite &pol1, APolarSite &pol2);
        inline vec  FieldPermESF(vec r, APolarSite &pol);
        inline void FieldIndu(APolarSite &pol1, APolarSite &pol2);
        inline void FieldInduAlpha(APolarSite &pol1, APolarSite &pol2);

        inline double PotentialPerm(vec r, APolarSite &pol);
        
        inline double E_f(APolarSite &pol1, APolarSite &pol2);
        inline double E_m(APolarSite &pol1, APolarSite &pol2);

        void ResetEnergy() { EP = EU_INTER = EU_INTRA = 0.0; 
                             EPP = EPU = EUU = 0.0;  }
        double &getEP() { return EP; }
        double &getEU_INTER() { return EU_INTER; }
        double &getEU_INTRA() { return EU_INTRA; }
        
        double &getEPP() { return EPP; }
        double &getEPU() { return EPU; }
        double &getEUU() { return EUU; }

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
        QMultipole     *_em;

    };


    // ++++++++++++++++++++++++++++++++++++++ //
    // Site workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class SiteOpQMultipole : public Thread
    {
    public:

        SiteOpQMultipole(int id, Topology *top,
                     QMultipole *master)
                   : _id(id), _top(top), _seg(NULL),
                     _master(master)
                   { _actor = Interactor(top,_master); };

       ~SiteOpQMultipole();

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void InitSlotData(Topology *top);
        void Run(void);

        void   EvalSite(Topology *top, Segment *seg);
        void   Charge(int state);
        int    Induce(int state);
        double Energy(int state);
        double EnergyStatic(int state);
        void   Depolarize();


    public:

        int                           _id;
        Topology                     *_top;
        Segment                      *_seg;
        QMultipole                  *_master;

        vector< Segment* >           _segsPolSphere; // Segments    in c/o 0-1
        vector< Segment* >           _segsOutSphere; // Segments    in c/0 1-2
        vector< vector<APolarSite*> > _polsPolSphere; // Polar sites in c/o 0-1
        vector< vector<APolarSite*> > _polsOutSphere; // Polar sites in c/o 1-2
        vector< vector<APolarSite*> > _APolarSites;    // Copy of top polar sites
        Interactor                   _actor;
    };


private:

    // Allocation of polar sites to fragments and segments
    map<string, vector<APolarSite*> >     _map_seg_APolarSites;
    map<string, vector<bool> >           _map_seg_chrgStates;
    map<string, vector<int> >            _alloc_frag_mpoleIdx;
    map<string, vector<string> >         _alloc_frag_mpoleName;
    map<string, bool>                    _map2md;

    // Thread management
    vector<Segment*> ::iterator _nextSite;
    vector<Segment*>            _listSegs;
    Mutex                       _nextSiteMutex;
    Mutex                       _coutMutex;
    Mutex                       _logMutex;
    bool                        _maverick;

    // Control options
    bool            _induce;
    int             _firstSeg;
    int             _lastSeg;
    vector<int>     _listSegIds;
    
    string          _checkPolesPDB;

    // ESP
    bool            _calcESP;
    bool            _ESPdoSystem;
    string          _espCubeFile;
    string          _espOutFile;
    vector<vec>     _cube;

    // ESF
    bool            _calcESF;
    string          _esfGridFile;
    string          _esfOutFile;
    vector<vec>     _esfGrid;

    // MolPol
    bool            _calcAlphaMol;
    bool            _doSystem;
    string          _alphaOutFile;

    // Logging
    string          _outFile;
    bool            _energies2File;
    map<int,vector<int> > _log_seg_iter;
    map<int,int>    _log_seg_sphereSize;
    map<int,vec>    _log_seg_com;
    vector<Segment*> _log_evaluated_segs;

    // Interaction parameters
    bool            _useCutoff;
    double          _cutoff;
    double          _cutoff2;
    bool            _useExp;
    double          _aDamp;
    bool            _useScaling;
    vector<double>  _scale1;

    // Convergence parameters
    float           _wSOR_N;
    float           _wSOR_C;
    double          _epsTol;
    int             _maxIter;

};


/**
 * Loads options and parameters from XML File including:
 * ... Thole model parameters
 * ... SOR parameters (convergence)
 * ... Control options (first, last seg., ...)
 */
void QMultipole::Initialize(Property *opt) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;

    cout << endl <<  "... ... Parametrizing Thole model";

    string key;
    string xmlfile;

    /* ---- OPTIONS.XML Structure ----
     * <qmultipole>
     *
     *      <multipoles></multipoles>
     *
     *      <control>
     *          <induce></induce>
     *          <first></first>
     *          <last></last>
     *          <output></output>
     *          <check></check>
     *      </control>
     *
     *      <esp>
     *          <calcESP></calcESP>
     *          <doSystem></doSystem>
     *          <cube></cube>
     *          <output></output>
     *      <esp>
     *
     *      <esf>
     *          <calcESF></calcESF>
     *          <grid></grid>
     *          <output></output>
     *      <esf>
     *
     *      <alphamol>
     *          <calcAlpha></calcAlpha>
     *          <doSystem></doSystem>
     *          <output></output>
     *      </alphamol>
     *
     *      <tholeparam>
     *          <cutoff></cutoff>
     *          <cutoff2></cutoff2>
     *          <expdamp></expdamp>
     *          <scaling></scaling>
     *      </tholeparam>
     *
     *      <convparam>
     *          <wSOR_N></wSOR_N>
     *          <wSOR_C></wSOR_C>
     *          <maxiter></maxiter>
     *          <tolerance></tolerance>
     *          <log></log>
     *      </convparam>
     */

    key = "options.qmultipole.multipoles";

        if ( opt->exists(key) ) {
            xmlfile = opt->get(key).as< string >();
        }

    key = "options.qmultipole.control";

        if ( opt->exists(key+".induce") ) {
            int induce = opt->get(key+".induce").as< int >();
            _induce = (induce == 0) ? false : true;
        }
        else { _induce = true; }

        if ( opt->exists(key+".first") ) {
            _firstSeg = opt->get(key+".first").as< int >();
        }
        else { _firstSeg = 1; }

        if ( opt->exists(key+".last") ) {
            _lastSeg = opt->get(key+".last").as< int >();
        }
        else { _lastSeg = -1; }

        if ( opt->exists(key+".output") ) {
            _outFile = opt->get(key+".output").as< string >();
            _energies2File = true;
        }
        else { _energies2File = false; }

        if ( opt->exists(key+".check") ) {
            _checkPolesPDB = opt->get(key+".check").as< string >();
        }
        else { _checkPolesPDB = ""; }

        if ( opt->exists(key+".list") ) {
            _listSegIds = opt->get(key+".list").as< vector<int> >();
        }
        else { ; }

    key = "options.qmultipole.esp";

        if ( opt->exists(key+".calcESP") ) {
            int calcESP = opt->get(key+".calcESP").as< int >();
            _calcESP = (calcESP == 0) ? false : true;
        }
        else {
            _calcESP = false;
        }
        if (_calcESP) {
            _espCubeFile = opt->get(key+".cube").as< string >();
            _espOutFile = opt->get(key+".output").as< string >();
            if (opt->exists(key+".doSystem")) {
                int doSystem = opt->get(key+".doSystem").as< int >();
                _ESPdoSystem = (doSystem == 0) ? false : true;
            }
            else { _ESPdoSystem = false; }
        }

    key = "options.qmultipole.esf";

        if ( opt->exists(key+".calcESF") ) {
            int calcESF = opt->get(key+".calcESF").as< int >();
            _calcESF = (calcESF == 0) ? false : true;
        }
        else {
            _calcESF = false;
        }
        if (_calcESF) {
            _esfGridFile = opt->get(key+".grid").as< string >();
            _esfOutFile = opt->get(key+".output").as< string >();
        }

    key = "options.qmultipole.alphamol";

        if ( opt->exists(key+".calcAlpha") ) {
            int calcAlpha = opt->get(key+".calcAlpha").as< int >();
            _calcAlphaMol = (calcAlpha == 0) ? false : true;
        }
        if ( opt->exists(key+".doSystem") ) {
            int doSystem = opt->get(key+".doSystem").as< int >();
            _doSystem = (doSystem == 0) ? false : true;
        }
        else {
            _doSystem = false;
        }
        if ( opt->exists(key+".output") ) {
            string alphaOutFile = opt->get(key+".output").as< string >();
            _alphaOutFile = alphaOutFile;
        }
        else {
            _alphaOutFile = "nofile";
        }

    key = "options.qmultipole.tholeparam";

        if ( opt->exists(key+".cutoff") ) {
            _cutoff = opt->get(key+".cutoff").as< double >();
            if (_cutoff) { _useCutoff = true; }
        }
        if ( opt->exists(key+".cutoff2") ) {
            _cutoff2 = opt->get(key+".cutoff2").as< double >();
        }
        else {
            _cutoff2 = _cutoff;
        }
        if ( opt->exists(key+".expdamp") ) {
            _aDamp = opt->get(key+".expdamp").as< double >();
            if (_aDamp) { _useExp = true; }
        }
         if ( opt->exists(key+".scaling") ) {
            _scale1 = opt->get(key+".scaling").as< vector<double> >();
            if (0 < _scale1.size() && _scale1.size() < 4) {
                _useScaling = true; }
            else {
                _useScaling = false;
                cout << endl << "... ... WARNING: 1-N SCALING SWITCHED OFF"; }
        }

    key = "options.qmultipole.convparam";

        if ( opt->exists(key+".wSOR_N") ) {
            _wSOR_N = opt->get(key+".wSOR_N").as< float >();
        }
        else { _wSOR_N = 0.75; }
        if ( opt->exists(key+".wSOR_C") ) {
            _wSOR_C = opt->get(key+".wSOR_C").as< float >();
        }
        else { _wSOR_C = 0.75; }

        if ( opt->exists(key+".maxiter") ) {
            _maxIter = opt->get(key+".maxiter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; }

    //if (!top->isEStatified()) { this->EStatify(top, opt); }

    if (_calcESP && (!_ESPdoSystem))      { this->CalculateESPInput(NULL); }
    if (this->_calcESF)                   { this->CalculateESF(NULL); }
    if (this->_calcAlphaMol)              { this->CalculateAlphaInput(NULL); }
}


/**
 * Creates polar-site templates using mapping file
 * ... Establishes allocation of polar sites to fragments
 * ... Calls GDMA punch-file parser to load multipoles and polarizabilities
 * ... In doing this, establishes which charge states are available (N, C, A)
 */
void QMultipole::EStatify(Topology *top, Property *options) {

    cout << endl << "... ... Estatify system: ";

    string key = "options.qmultipole";
    string allocFile = options->get(key+".multipoles").as<string> ();

    // ++++++++++++++++++++++++++++++++ //
    // Load polar-site indices from XML //
    // ++++++++++++++++++++++++++++++++ //

    // => Output to maps:
    map<string, vector<int> > alloc_frag_mpoleIdx;
    map<string, vector<string> > alloc_frag_mpoleName;
    map<string, string > alloc_seg_mpoleFiles_n;
    map<string, string > alloc_seg_mpoleFiles_e;
    map<string, string > alloc_seg_mpoleFiles_h;

    Property allocation; // <- Which polar sites are part of which fragment?
    load_property_from_xml(allocation, allocFile.c_str());


    /* --- MULTIPOLES.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name>DCV</name>
     *
     *              <multipoles_n></multipoles_n>
     *              <multipoles_e></multipoles_e>
     *              <multipoles_h></multipoles_h>
     *
     *              <map2md></map2md>
     *
     *              <fragments>
     *                  <fragment>
     *                  <name></name>
     *                  <mpoles></mpoles>
     *                  </fragment>
     *              </fragments>
     *              ...
     *              ...
     */


    key = "topology.molecules.molecule";
    list<Property *> mols = allocation.Select(key);
    list<Property *>::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); molit++) {

        string molName = (*molit)->get("name").as<string> ();

        key = "segments.segment";
        list<Property *> segs = (*molit)->Select(key);
        list<Property *>::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); segit++) {

            string segName = (*segit)->get("name").as<string> ();

            // GDMA filenames for cation (h), neutral (n), anion (e) state
            string mpoleFile_n = (*segit)->get("multipoles_n").as<string> ();
            alloc_seg_mpoleFiles_n[segName] = mpoleFile_n;
            if ( (*segit)->exists("multipoles_e")) {
                string mpoleFile_e = (*segit)->get("multipoles_e").as<string>();
                alloc_seg_mpoleFiles_e[segName] = mpoleFile_e;
            }
            if ( (*segit)->exists("multipoles_h")) {
                string mpoleFile_h = (*segit)->get("multipoles_h").as<string>();
                alloc_seg_mpoleFiles_h[segName] = mpoleFile_h;
            }

            // Default: Project multipoles onto rigidified coordinates
            if ( (*segit)->exists("map2md")) {
                int map2md = (*segit)->get("map2md").as<int>();
                _map2md[segName] = (map2md == 0) ? false : true;
            }
            else {
                _map2md[segName] = false; // i.e. map to rigidified coordinates
            }

            key = "fragments.fragment";
            list<Property *> frags = (*segit)->Select(key);
            list<Property *>::iterator fragit;

            for (fragit = frags.begin(); fragit != frags.end(); fragit++) {

                string fragName = (*fragit)->get("name").as<string> ();
                string mapKeyName = fragName + segName + molName;

                string mpoles = (*fragit)->get("mpoles").as<string> ();

                Tokenizer tokPoles(mpoles, " \t\n");
                vector<string> mpoleInfo;
                tokPoles.ToVector(mpoleInfo);

                vector<int> mpoleIdcs;
                vector<string> mpoleNames;

                vector<string> ::iterator strit;
                for (strit=mpoleInfo.begin(); strit<mpoleInfo.end(); strit++) {

                    Tokenizer tokPoleInfo( (*strit), " :");
                    vector<string> poleInfo;
                    tokPoleInfo.ToVector(poleInfo);

                    int mpoleIdx = boost::lexical_cast<int>(poleInfo[0]);
                    string mpoleName = poleInfo[1];

                    mpoleIdcs.push_back(mpoleIdx);
                    mpoleNames.push_back(mpoleName);

                }

                alloc_frag_mpoleIdx[mapKeyName] = mpoleIdcs;
                alloc_frag_mpoleName[mapKeyName] = mpoleNames;
            }
        }
    }

    // ++++++++++++++++++++++ //
    // Parse GDMA punch files //
    // ++++++++++++++++++++++ //

    // => Output to APolarSite Container (template container)
    map<string, vector<APolarSite*> > map_seg_APolarSites;
    map<string, vector<bool> >       map_seg_chrgStates;


    // Multipoles for neutral state
    int state = 0;
    map<string, string > ::iterator strit;
    for (strit = alloc_seg_mpoleFiles_n.begin();
         strit != alloc_seg_mpoleFiles_n.end();
         strit++) {

        string segName = strit->first;
        string filename = strit->second;

        vector< APolarSite* > poles = ParseGdmaFile(filename, state);

        map_seg_APolarSites[segName] = poles;
        map_seg_chrgStates[segName].resize(3);
        map_seg_chrgStates[segName][0] = false; // <- negative
        map_seg_chrgStates[segName][1] = true;  // <- neutral
        map_seg_chrgStates[segName][2] = false; // <- positive
        poles.clear();

    }

    // Multipoles for negative state
    state = -1;
    if ( alloc_seg_mpoleFiles_e.size() ) {
        for (strit = alloc_seg_mpoleFiles_e.begin();
             strit != alloc_seg_mpoleFiles_e.end();
             strit++) {
            string segName = strit->first;
            string filename = strit->second;

            vector< APolarSite* > polesAnion = ParseGdmaFile(filename, state);
            map_seg_chrgStates[segName][state+1] = true;

            // Merge with polar sites for neutral state
            vector< APolarSite* > polesNeutral = map_seg_APolarSites[segName];

            assert(polesAnion.size() == polesNeutral.size());
            for (unsigned int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesAnion[i]->getQs(state), state );
                polesNeutral[i]->setPs( polesAnion[i]->getPs(state), state );
                delete polesAnion[i];
            }
            polesAnion.clear();
        }
    }

    // Multipoles for positive state
    state = +1;
    if ( alloc_seg_mpoleFiles_h.size() ) {
        for (strit = alloc_seg_mpoleFiles_h.begin();
             strit != alloc_seg_mpoleFiles_h.end();
             strit++) {
            string segName = strit->first;
            string filename = strit->second;

            vector< APolarSite* > polesCation = ParseGdmaFile(filename, state);
            map_seg_chrgStates[segName][state+1] = true;

            // Merge with polar sites for neutral state
            vector< APolarSite* > polesNeutral = map_seg_APolarSites[segName];

            assert(polesCation.size() == polesNeutral.size());
            for (unsigned int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesCation[i]->getQs(state), state );
                polesNeutral[i]->setPs( polesCation[i]->getPs(state), state );
                delete polesCation[i];
            }
            polesCation.clear();
        }
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Forward information on polar-site templates to topology //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    // Containers to use for mapping
    // map<string, vector<APolarSite*> >     map_seg_APolarSites;
    // map<string, vector<int> >            alloc_frag_mpoleIdx;
    // map<string, vector<string> >         alloc_frag_mpoleName;

    // TO CONTINUE FROM HERE ...
    // Loop over segments, loop over fragments in segments.
    // => For each segment, look up associated polar sites using map 1;
    //      => For each fragment, look up associated
    //         polar sites idcs. using map 2;
    //          => For each index in retrieved container, create new polar site,
    //             importing information from polar site container of segment
    //          => Translate positions to MD frame, copy local frame from
    //             fragment

    _map_seg_APolarSites = map_seg_APolarSites;
    _map_seg_chrgStates = map_seg_chrgStates;
    _alloc_frag_mpoleIdx =  alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName = alloc_frag_mpoleName;

}


/**
 * Parses GDMA punch file (for structure sample see below)
 * ... Loads positions, rank, multipoles
 * ... Converts from a.u. to internal units
 * ... Determines polarizabilities from element type
 * @return Yields polar-site template container
 */
vector<APolarSite*> QMultipole::ParseGdmaFile(string filename, int state) {

    return APS_FROM_MPS(filename,state);    
    
}


/**
 * Parses GAUSSIAN cube file (for structure sample see below)
 * ... NOTE Only parses header to extract voxel vectors and repetition numbers
 */
vector<vec> QMultipole::ParseCubeFileHeader(string filename) {

    vector< vec > cube;
    cube.resize(5);

    //               STRUCTURE CUBE CONTAINER
    //       0        1         2         3         4
    //  |         |         |         |         |         |
    //   ox oy oz  Nx Ny Nz  xx xy xz  yx yy yz  zx zy zz

    std::string line;
    int lineCount = 0;

    double Nx = 0;
    double Ny = 0;
    double Nz = 0;
    double BOHR2NM = 0.0529189379;
    double ANGSTROM2NM = 0.1;
    double ext2nm = BOHR2NM;

    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open()) {
        while (intt.good()) {

            ++lineCount;

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            // First two lines are comment lines
            if (lineCount == 1 || lineCount == 2) { continue; }



            // Skip empty lines
            if ( split.size() == 0) { continue; }

            // => Interesting information here
            if (lineCount == 3) {

                assert( split.size() == 4 );
                double ox = boost::lexical_cast<double>( split[1] );
                double oy = boost::lexical_cast<double>( split[2] );
                double oz = boost::lexical_cast<double>( split[3] );

                cube[0] = vec(ox, oy, oz);
            }

            /*
            STRUCTURE SAMPLE CPMD CUBE FILE
            OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
              3    0.000000    0.000000    0.000000
             40    0.283459    0.000000    0.000000
             40    0.000000    0.283459    0.000000
             40    0.000000    0.000000    0.283459
              8    0.000000    5.570575    5.669178    5.593517
              1    0.000000    5.562867    5.669178    7.428055
              1    0.000000    7.340606    5.669178    5.111259
              :             :          :           :            :            :
              :             :          :           :            :            :
              :             :          :           :            :            :
                In this case there will be 40 x 40 x 40 floating point values
              :             :          :           :            :            :
              :             :          :           :            :            :
              :             :          :           :            :            :
            */

            if (lineCount == 4) {

                assert( split.size() == 4 );
                Nx = boost::lexical_cast<double>( split[0] );
                ext2nm = (Nx > 0) ? BOHR2NM : ANGSTROM2NM;
                Nx = (Nx >= 0) ? Nx : -Nx;
                double xx = ext2nm * boost::lexical_cast<double>( split[1] );
                double xy = ext2nm * boost::lexical_cast<double>( split[2] );
                double xz = ext2nm * boost::lexical_cast<double>( split[3] );

                cube[2] = vec(xx, xy, xz);
            }

            if (lineCount == 5) {

                assert( split.size() == 4 );
                Ny = boost::lexical_cast<double>( split[0] );
                ext2nm = (Ny > 0) ? BOHR2NM : ANGSTROM2NM;
                Ny = (Ny >= 0) ? Ny : -Ny;
                double yx = ext2nm * boost::lexical_cast<double>( split[1] );
                double yy = ext2nm * boost::lexical_cast<double>( split[2] );
                double yz = ext2nm * boost::lexical_cast<double>( split[3] );

                cube[3] = vec(yx, yy, yz);
            }

            if (lineCount == 6) {

                assert( split.size() == 4 );
                Nz = boost::lexical_cast<double>( split[0] );
                ext2nm = (Nz > 0) ? BOHR2NM : ANGSTROM2NM;
                Nz = (Nz >= 0) ? Nz : -Nz;
                double zx = ext2nm * boost::lexical_cast<double>( split[1] );
                double zy = ext2nm * boost::lexical_cast<double>( split[2] );
                double zz = ext2nm * boost::lexical_cast<double>( split[3] );

                cube[4] = vec(zx, zy, zz);
                cube[1] = vec(Nx, Ny, Nz);
            }

            // Header at finish?
            if (lineCount > 6) {
                break;
            }
        }
    }
    else { cout << endl << "ERROR: No such file " << filename << endl; }

    // Convert origin
    cube[0] *= ext2nm;

    return cube;
}


/**
 * Calculates electrostatic potential at selected grid points
 * ... NOTE Calculation done for input coordinates (as found in GDMA file)
 * ... NOTE Grid is read in from file in Gaussian cube format
 */
void QMultipole::CalculateESPInput(Topology *top) {

    cout << endl << "... ... Calculating ESP";

    // Load grid from cube file
    _cube = ParseCubeFileHeader(_espCubeFile);

    // +++++++++++++++++++++++++++++++++++++++++++ //
    // Iterate over segment types (in GDMA config) //
    // +++++++++++++++++++++++++++++++++++++++++++ //

    map< string, vector< APolarSite* > >::iterator sit;
    for (sit = this->_map_seg_APolarSites.begin();
         sit != this->_map_seg_APolarSites.end();
         ++sit) {

        cout << endl << "... ... ... ESP for segment type " << sit->first
             << flush;

        vector< APolarSite* > poles = sit->second;
        vector< APolarSite* >::iterator pit;

        // Translate CoM to origin
        vec CoM = vec(0.,0.,0.);
        int N = poles.size();
        for (pit = poles.begin(); pit < poles.end(); ++pit) {
            CoM += (*pit)->getPos();
        }
        CoM /= N;
        for (pit = poles.begin(); pit < poles.end(); ++pit) {
            vec newpos = (*pit)->getPos() - CoM;
            (*pit)->setPos( newpos );
        }

        for (int state = -1; state < 2; ++state) {

            // Multipole data for this state?
            bool hasState = _map_seg_chrgStates[sit->first][state+1];
            if (!hasState) { continue; }

            // Charge container appropriately
            for (pit = poles.begin(); pit < poles.end(); ++pit) {
                (*pit)->Depolarize();
                (*pit)->Charge(state);
            }

            // Prepare cube file, for structure sample see cube-file parser
            FILE *out;
            string espOutFile =  sit->first
                            + "_" + boost::lexical_cast<string>(state)
                            + "_" + _espOutFile;
            out = fopen(espOutFile.c_str(), "w");

            fprintf(out, "CUBE FILE GENERATED BY VOTCA XTP \n");
            fprintf(out, "ESP %4s CHARGE %1d \n", (sit->first).c_str(), state);

            cout << endl << "... ... ... ... State " << state << ": " << flush;

            this->CalculateESP(poles, out);

            fclose(out);
        }
    }
}


/**
 * Calculates electrostatic potential at selected grid points
 * ... NOTE Calculation done for rigid coordinates (as found in trajectory)
 * ... NOTE Grid is read in from file in Gaussian cube format
 */
void QMultipole::CalculateESPRigid(Topology *top) {

    cout << endl << "... ... Calculate ESP for segment ID: " << flush;

    // Cube volumetric grid
    _cube = ParseCubeFileHeader(_espCubeFile);

    int segId;
    cin >> segId;

    Segment *seg = top->getSegment(segId);

    cout << "... ... ... Selected segment type " << seg->getName() << flush;

    vector< APolarSite* > poles = seg->APolarSites();
    vector< APolarSite* > ::iterator pit;

    // Translate CoM to origin
    vec CoM = vec(0.,0.,0.);
    int N = poles.size();
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        CoM += (*pit)->getPos();
    }
    CoM /= N;
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        vec newpos = (*pit)->getPos() - CoM;
        (*pit)->setPos( newpos );
    }

    // Investigate charge states
    for (int state = -1; state < 2; ++state) {
        bool hasState = _map_seg_chrgStates[seg->getName()][state+1];
        if (!hasState) { continue; }

        for (pit = poles.begin(); pit < poles.end(); ++pit) {
            (*pit)->Charge(state);
        }

        FILE *out;
        string espOutFile = seg->getName()
                        + "_ID" + boost::lexical_cast<string>(seg->getId())
                        + "_" + boost::lexical_cast<string>(state)
                        + "_" + _espOutFile;
        out = fopen(espOutFile.c_str(), "w");

        fprintf(out, "CUBE FILE GENERATED BY VOTCA XTP \n");
        fprintf(out, "ESP %4s ID %4d CHARGE %1d \n",
                      seg->getName().c_str(), seg->getId(), state);

        cout << endl << "... ... ... State " << state << ": " << flush;

        this->CalculateESP(poles, out);

        fclose(out);
    }
}


/**
 * Calculates ESP for arbitrary assembly of polar sites
 * ... NOTE Writes output in cube-file format, two header lines should
 * ...      already have been written to FILE* out
 */
void QMultipole::CalculateESP(vector< APolarSite* > &poles, FILE *out) {

    // Calculate potential for each voxel, write to output on the fly
    Interactor actor = Interactor(NULL, this);
    vector< APolarSite* >::iterator pit;

    vec Cxyz = _cube[0];

    vec C0 = _cube[0];
    vec Cx = _cube[2];
    vec Cy = _cube[3];
    vec Cz = _cube[4];

    int Nx = int(_cube[1].getX() + 0.5);
    int Ny = int(_cube[1].getY() + 0.5);
    int Nz = int(_cube[1].getZ() + 0.5);


    // Prepare cube file: Header
    double NM2BOHR = 1 / 0.0529189379;
    fprintf(out, "%5ld %3.7f %3.7f %3.7f \n",
                  poles.size(),
                  Cxyz.getX()*NM2BOHR,
                  Cxyz.getY()*NM2BOHR,
                  Cxyz.getZ()*NM2BOHR);
    fprintf(out, "%5d %3.7f %3.7f %3.7f \n",
                  Nx,
                  Cx.getX()*NM2BOHR,
                  Cx.getY()*NM2BOHR,
                  Cx.getZ()*NM2BOHR);
    fprintf(out, "%5d %3.7f %3.7f %3.7f \n",
                  Ny,
                  Cy.getX()*NM2BOHR,
                  Cy.getY()*NM2BOHR,
                  Cy.getZ()*NM2BOHR);
    fprintf(out, "%5d %3.7f %3.7f %3.7f \n",
                  Nz,
                  Cz.getX()*NM2BOHR,
                  Cz.getY()*NM2BOHR,
                  Cz.getZ()*NM2BOHR);
    for (pit = poles.begin(); pit < poles.end(); ++pit) {

        string elem = (*pit)->_name;
        int Z = 0;
        if      (elem == "C") { Z = 6; }
        else if (elem == "H") { Z = 1; }
        else if (elem == "N") { Z = 7; }
        else if (elem == "O") { Z = 8; }
        else if (elem == "S") { Z = 16; }

        fprintf(out, "%5d %3.7f %3.7f %3.7f %3.7f \n",
                Z,
                0.000,
                (*pit)->getPos().getX() * NM2BOHR,
                (*pit)->getPos().getY() * NM2BOHR,
                (*pit)->getPos().getZ() * NM2BOHR);
    }

    int phiCount = 0;
    double int2V = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for ( int k = 0; k < Nz; ++k) {

                Cxyz = C0 + i*Cx + j*Cy + k*Cz;

                double phi = 0.0;
                ++phiCount;

                for (pit = poles.begin(); pit < poles.end(); ++pit) {
                    phi += actor.PotentialPerm(Cxyz, *(*pit));
                }

                fprintf(out, " %1.7E", phi * int2V);
                if ( (phiCount % 6) == 0 ) {
                    fprintf(out, " \n");
                    phiCount = 0;
                }
            }
        }
    }

    cout << " Calculated potential at " << Nx*Ny*Nz
         << " grid points. " << flush;

    // Calculate self-energy of configuration (just out of interest)
    double E = 0.0;
    vector< APolarSite* > ::iterator pit2;
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        for (pit2 = pit + 1; pit2 < poles.end(); ++pit2) {
            E += actor.EnergyInterESP(*(*pit), *(*pit2));
        }
    }
    double E2eV = 1.602176487e-19 * 1e9 * 1/(4*M_PI*8.854187817e-12);
    cout << endl << "... ... ... ... "
         << "Electrostatic self-energy of configuration: "
         << E*E2eV << " eV. ";

}


/**
 * Calculates electrostatic field at selected grid points
 * ... NOTE Calculation done for rigid coordinates (as found in GDMA file)
 * ... NOTE Grid is currently loaded from xyz file
 */
void QMultipole::CalculateESF(Topology *top) {

    cout << endl << "... ... Calculating ESF";

    //double int2N_C = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 * 1.000e-18;

    // +++++++++++++++++++ //
    // Load grid from file //
    // +++++++++++++++++++ //

    vector<vec> gridPoints;
    double CONVERT2NM = 1;

    std::string line;
    std::ifstream intt;
    intt.open(this->_esfGridFile.c_str());

    if (intt.is_open() ) {

        while ( intt.good() ) {


            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }

            // ! Interesting information here, e.g.
            // Units bohr

            //  -0.3853409355   0.0004894995  -0.0003833545
            //  -0.0002321905   0.2401559510   0.6602334308
            // ! ...

            // Units used
            if ( split[0] == "Units") {
                string units = split[1];

                if      (units == "bohr") { CONVERT2NM = 0.0529189379; }
                else if (units == "nm") { CONVERT2NM = 1.; }
                else if (units == "A") { CONVERT2NM = 0.1; }
                else {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + _esfGridFile + " not supported.");
                }
            }

            else {

                double x = CONVERT2NM * boost::lexical_cast<double>(split[0]);
                double y = CONVERT2NM * boost::lexical_cast<double>(split[1]);
                double z = CONVERT2NM * boost::lexical_cast<double>(split[2]);

                gridPoints.push_back(vec(x,y,z));
            }
        }
    }
    else { cout << endl << "ERROR: No such file " << _esfGridFile << endl; }

    this->_esfGrid = gridPoints;

    // +++++++++++++++++++++++++++++++++++++++++++ //
    // Iterate over segment types (in GDMA config) //
    // +++++++++++++++++++++++++++++++++++++++++++ //

    map< string, vector< APolarSite* > >::iterator sit;
    for (sit = this->_map_seg_APolarSites.begin();
         sit != this->_map_seg_APolarSites.end();
         ++sit) {

        cout << endl << "... ... ... ESF for segment " << sit->first << flush;



        vector< APolarSite* > poles = sit->second;
        vector< APolarSite* >::iterator pit;

        for (int state = -1; state < 2; ++state) {

            // Multipole data for this state?
            bool hasState = _map_seg_chrgStates[sit->first][state+1];
            if (!hasState) { continue; }

            // Store the field
            vector< vec > gridPointESF;

            // Charge container appropriately
            for (pit = poles.begin(); pit < poles.end(); ++pit) {
                (*pit)->Charge(state);
            }

            // Calculate field for each grid point, store in gridPointESF
            Interactor actor = Interactor(top, this);

            vector< vec >::iterator grit;
            for (grit = _esfGrid.begin(); grit < _esfGrid.end(); ++grit) {

                vec F = vec(0.0, 0.0, 0.0);

                for (pit = poles.begin(); pit < poles.end(); ++pit) {
                    F += actor.FieldPermESF((*grit), *(*pit));
                }
                gridPointESF.push_back(F);
            }

            cout << endl << "... ... ... ... State " << state << ": "
                         << " Calculated field at " << gridPointESF.size()
                         << " grid points. " << flush;


            // Write to file
            FILE *out;
            string esfOutFile =  sit->first
                            + "_" + boost::lexical_cast<string>(state)
                            + "_" + _esfOutFile;

            out = fopen(esfOutFile.c_str(), "w");

            assert(_esfGrid.size() == gridPointESF.size());
            unsigned int g = 0;
            unsigned int p = 0;
            for ( ; p < gridPointESF.size(); ++g, ++p) {

                fprintf(out, " %3.8f %3.8f %3.8f   %3.8f %3.8f %3.8f \n",
                        _esfGrid[g].getX(),
                        _esfGrid[g].getY(),
                        _esfGrid[g].getZ(),
                        gridPointESF[p].getX(),
                        gridPointESF[p].getY(),
                        gridPointESF[p].getZ() );
            }
        }
    }
}


/**
 * Calculates molecular polarizability tensor for GDMA coordinates
 * ... Uses Thole SCF model with exponential damping function
 */
void QMultipole::CalculateAlphaInput(Topology *top) {

    cout << endl << "... ... Calculating site polarizabilities";

    map< string, vector< APolarSite* > >::iterator sit;
    for (sit = _map_seg_APolarSites.begin();
         sit != _map_seg_APolarSites.end();
         ++sit) {

        cout << endl << "... ... ... Polarizability for segment "
             << sit->first << flush;

        FILE *out;
        string alphaOutFile = sit->first + "_" + _alphaOutFile;
        out = fopen(alphaOutFile.c_str(), "w");

        vector< APolarSite* > ::iterator pit;

        // Neutral state
        int state = 0;
        if (_map_seg_chrgStates[sit->first][state+1]) {

            for (pit = sit->second.begin(); pit < sit->second.end(); ++pit) {
                (*pit)->Charge(state);
            }

            cout << endl << "... ... ... ... NEUTRAL " << state;
            fprintf(out, "[ NEUTRAL ] \n");
            this->CalculateAlphaRigid(sit->second, out, false);
        }

        // Cationic state
        state = +1;
        if (_map_seg_chrgStates[sit->first][state+1]) {

            for (pit = sit->second.begin(); pit < sit->second.end(); ++pit) {
                (*pit)->Charge(state);
            }

            cout << endl << "... ... ... ... CATION " << state;
            fprintf(out, "\n\n[ CATION ] \n");
            this->CalculateAlphaRigid(sit->second, out, false);
        }

        // Anionic state
        state = -1;
        if (_map_seg_chrgStates[sit->first][state+1]) {

            for (pit = sit->second.begin(); pit < sit->second.end(); ++pit) {
                (*pit)->Charge(state);
            }

            cout << endl << "... ... ... ... ANION " << state;
            fprintf(out, "\n\n[ ANION ] \n");
            this->CalculateAlphaRigid(sit->second, out, false);
        }

        // Not necessary, still: Charge back to neutral state
        for (pit = sit->second.begin(); pit < sit->second.end(); ++pit) {
                (*pit)->Charge(0);
        }

    }
}


/**
 * Calculates polarizability for arbitrary assembly of polar sites
 */
void QMultipole::CalculateAlphaRigid(vector< APolarSite* > &poles,
                                      FILE *out, bool system) {

    vector< APolarSite* > ::iterator pit;

    double axx, axy, axz;             // |
    double ayx, ayy, ayz;             // |-> Polarizability tensor
    double azx, azy, azz;             // |

    double siteU1x, siteU1y, siteU1z; // |-> Molecular ind. dipole

    double extFx, extFy, extFz;       // |-> External applied field

    // +++++++++++++++++++++++++++++ //
    // External field in x-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.1;
    extFy = 0.0;
    extFz = 0.0;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    this->SiteAlphaInduce(NULL, poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
    }

    // Calculate associated column of polarizability tensor
    axx = - siteU1x / extFx;
    ayx = - siteU1y / extFx;
    azx = - siteU1z / extFx;

    // +++++++++++++++++++++++++++++ //
    // External field in y-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.0;
    extFy = 0.1;
    extFz = 0.0;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    this->SiteAlphaInduce(NULL, poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
    }

    // Calculate associated column of polarizability tensor
    axy = - siteU1x / extFy;
    ayy = - siteU1y / extFy;
    azy = - siteU1z / extFy;

    // +++++++++++++++++++++++++++++ //
    // External field in z-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.0;
    extFy = 0.0;
    extFz = 0.1;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    this->SiteAlphaInduce(NULL, poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
    }

    // Calculate associated column of polarizability tensor
    axz = - siteU1x / extFz;
    ayz = - siteU1y / extFz;
    azz = - siteU1z / extFz;

    // +++++++++++++++++++++++++++ //
    // Sum, Trace, Diagonalization //
    // +++++++++++++++++++++++++++ //

    // Sum over atomic polarizabilities
    double SUM_alpha = 0.0;
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        SUM_alpha += (*pit)->getIsoP();
    }

    double NM3_2_A3 = 1000.;
    SUM_alpha *= NM3_2_A3;

    // TODO Diagonalize polarizability tensor
    // For now, trace is enough
    double ISO_alpha = (axx + ayy + azz) / 3.;

    ISO_alpha *= NM3_2_A3;

    // Eigenvalues of polarizability tensor
    matrix alpha = matrix( vec(axx,ayx,azx),
                           vec(axy,ayy,azy),
                           vec(axz,ayz,azz) );
    matrix::eigensystem_t EIGEN;
    alpha.SolveEigensystem(EIGEN);


    if (! system) {
        cout << endl << "... ... ... ... Sum over atomic alphas:    "
             << SUM_alpha << " A**3.";
        cout << endl << "... ... ... ... 1/3 trace of alpha tensor: "
             << ISO_alpha << " A**3.";
        cout << endl << "... ... ... ... Eigenvalues: "
             << EIGEN.eigenvalues[0]*NM3_2_A3 << " A**3, "
             << EIGEN.eigenvalues[1]*NM3_2_A3 << " A**3, "
             << EIGEN.eigenvalues[2]*NM3_2_A3 << " A**3. ";
    }

    if (_alphaOutFile != "noname" && _alphaOutFile != "") {

        vec x = EIGEN.eigenvecs[0];
        vec y = EIGEN.eigenvecs[1];
        vec z = EIGEN.eigenvecs[2];

        if (! system) {
            fprintf(out, "Polarizability tensor in global frame \n");
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3);
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", ayx*NM3_2_A3, ayy*NM3_2_A3, ayz*NM3_2_A3);
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", azx*NM3_2_A3, azy*NM3_2_A3, azz*NM3_2_A3);
            fprintf(out, "\n");
            fprintf(out, "Polarizability tensor in eigenframe \n");
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", EIGEN.eigenvalues[0]*NM3_2_A3,0.,0.);
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", 0.,EIGEN.eigenvalues[1]*NM3_2_A3,0.);
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", 0.,0.,EIGEN.eigenvalues[2]*NM3_2_A3);
            fprintf(out, "\n");
            fprintf(out, "Polarizability tensor main axes \n");
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", x.getX(),x.getY(),x.getZ());
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", y.getX(),y.getY(),y.getZ());
            fprintf(out, "%4.7f  %4.7f  %4.7f \n", z.getX(),z.getY(),z.getZ());
        }
        else {
            fprintf(out, "%4.7f  %4.7f  %4.7f  %4.7f  %4.7f ",
                         ISO_alpha, SUM_alpha,
                         EIGEN.eigenvalues[0]*NM3_2_A3,
                         EIGEN.eigenvalues[1]*NM3_2_A3,
                         EIGEN.eigenvalues[2]*NM3_2_A3);
        }
    }
}


/**
 * Induces dipoles on polar sites in self-consistent manner until converged
 * ... NOTE Designed to be used in conjunction with ::CalculateAlphaMol()
 */
void QMultipole::SiteAlphaInduce(Topology *top, vector<APolarSite*> &poles) {

    int    maxIter = this->_maxIter;
    double wSOR = this->_wSOR_N;
    double eTOL = this->_epsTol;

    vector< APolarSite* >::iterator pit1;
    vector< APolarSite* >::iterator pit2;

    Interactor actor = Interactor(top, this);

    // Induce to first order
    for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
        (*pit1)->InduceDirect();
    }


    int iter = 0;
    for ( ; iter < maxIter; ++iter) {

        // Reset induction field
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            (*pit1)->ResetFieldU();
        }

        // Calculate higher-order induction field
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            for (pit2 = pit1 + 1; pit2 < poles.end(); ++pit2) {
                actor.FieldInduAlpha(*(*pit1),*(*pit2));
            }
        }

        // Induce dipoles
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            (*pit1)->Induce(wSOR);
        }

        // Check for convergence
        bool converged = true;
        double maxdU = -1;
        double avgdU = 0.0;
        int    baseN = 0;

        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
             double dU = (*pit1)->HistdU();
             avgdU += dU;
             ++baseN;
             if ( dU > maxdU ) { maxdU = dU; }
             if ( dU > eTOL ) { converged = false; }
        }

        avgdU /= baseN;
        if (avgdU < eTOL/10.) { converged = true; }

        //cout << endl << "ITER" << iter << " | MAX dU " << maxdU
        //     << " | AVG dU " << avgdU;

        // Break if converged
        if (converged) { break; }
        else if (iter == maxIter - 1) {
            cout << endl << "... ... ... "
                 << "WARNING Induced multipoles did not converge to precision: "
                 << " AVG dU:U " << avgdU << flush;
            break;
        }
    }
}


/**
 * Equips segments with polar sites using polar-site template container
 * ... Loops over segments + fragments in topology
 * ... Looks up associated polar sites, creates 'new' instances
 * ... Translates + rotates polar-site into global frame
 */
void QMultipole::DistributeMpoles(Topology *top) {

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg = *sit;
        vector<APolarSite*> poleSites = _map_seg_APolarSites.at(seg->getName());
        seg->setChrgStates(_map_seg_chrgStates[seg->getName()]);

        bool map2md = _map2md[seg->getName()];

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            string idkey = frag->getName() + seg->getName()
                         + seg->getMolecule()->getName();
            vector<int> polesInFrag = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag = _alloc_frag_mpoleName.at(idkey);

            if (map2md) {
                if (polesInFrag.size() != frag->Atoms().size()) {
                    cout << endl << "ERROR: Segment " << seg->getName() <<
                            " Fragment " << frag->getName() <<
                            ": MAP2MD = TRUE requires same number of polar sites "
                            " as there are atoms to perform mapping. " << endl;
                    throw runtime_error("User not paying attention. ");
                }
            }


            for (unsigned int i = 0; i < polesInFrag.size(); i++) {

                string name = namesInFrag[i];
                int poleId = polesInFrag[i];

                APolarSite *templ = poleSites[poleId-1];
                APolarSite *newSite = top->AddAPolarSite(name);
                newSite->ImportFrom(templ);
                seg->AddAPolarSite(newSite);
                frag->AddAPolarSite(newSite);

                // Shift + rotate
                if (!map2md) {
                    newSite->Translate(frag->getTransQM2MD());
                    newSite->Rotate(frag->getRotQM2MD(), frag->getCoMD());
                }
                else {
                    vec mdpos = frag->Atoms()[i]->getPos();
                    newSite->setPos(mdpos);
                    if (newSite->getRank() > 0) {
                        cout << endl << "ERROR: MAP2MD = TRUE "
                        " prevents use of higher-rank multipoles. " << endl;
                        throw runtime_error("User not paying attention. ");
                    }
                }
            }
        }
    }

    top->setIsEStatified(true);

}


/**
 * Creates threads to work on individual segments
 * ... Rigidifies + equips topology with polar sites
 * ... Creates, starts threads
 * ... Waits for threads to finish
 */
bool QMultipole::EvaluateFrame(Topology *top) {

    // ++++++++++++++++++++ //
    // Ridigidfy + Estatify //
    // ++++++++++++++++++++ //

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }

    // Forward multipoles to topology
    if (top->isEStatified() == false) {
        this->DistributeMpoles(top);
        cout << endl << "... ... Created " << top->APolarSites().size()
                 << " multipole sites." << flush;
    }
    else { cout << endl << "... ... System is already estatified."; }

    if (_checkPolesPDB != "") {
        cout << endl << "... ... Writing polar-site coordinates to "
             << _checkPolesPDB << ". " << flush;
        // To check rotation into global frame
        string mpNAME = _checkPolesPDB;
        FILE *mpPDB = NULL;
        mpPDB = fopen(mpNAME.c_str(), "w");
        vector<Segment*>::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {
            (*sit)->WritePDB(mpPDB, "Multipoles", "Charges");
        }
        fclose(mpPDB);

        string paNAME = "pa_check.xyz";
        FILE *paXYZ = NULL;
        paXYZ = fopen(paNAME.c_str(), "w");
        fprintf(paXYZ, "%6ld\n\n", 
                    4*top->getSegment(1)->APolarSites().size()
                   +4*top->getSegment(2)->APolarSites().size());
        vector<APolarSite*>::iterator pit;
        for (pit = top->getSegment(1)->APolarSites().begin();
             pit < top->getSegment(1)->APolarSites().end();
             ++pit) {
            (*pit)->PrintTensorPDB(paXYZ, -1);
        }
        for (pit = top->getSegment(2)->APolarSites().begin();
             pit < top->getSegment(2)->APolarSites().end();
             ++pit) {
            (*pit)->PrintTensorPDB(paXYZ, -1);
        }
        fclose(paXYZ);
    }

    if (_calcESP && _ESPdoSystem) {
        this->CalculateESPRigid(top);
    }

    // What has already been done? If anything at all, return 0.
    if ( (_calcAlphaMol && !_doSystem) || _calcESP || _calcESF) {

        cout << endl << "... ... Calculated ";

        if (this->_calcAlphaMol) {
            cout << "site polarizability tensors, ";
        }
        if (this->_calcESP) {
            cout << "ESPs, ";
        }
        if (this->_calcESF) {
            cout << "ESFs, ";
        }
        cout << "return. ";
        return 0;
    }

    // ++++++++++++++++++++++++++++++++++++++++ //
    // Polarizabilities for all sites in system //
    // ++++++++++++++++++++++++++++++++++++++++ //

    if (_calcAlphaMol && _doSystem) {

        cout << endl << "... ... Calculating site polarizabilities" << endl;

        FILE *out;
        string alphaOutFile =
                "Frame" + boost::lexical_cast<string>(top->getDatabaseId())
                + "_" + _alphaOutFile;
        out = fopen(alphaOutFile.c_str(), "w");

        vector< Segment* > ::iterator sit;
        vector< APolarSite* > ::iterator pit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {

            vector< APolarSite* > poles = (*sit)->APolarSites();


            cout << "\r... ... ... Segment " << (*sit)->getId() << " " << flush;

            fprintf(out, "%4d %4s ", (*sit)->getId(),(*sit)->getName().c_str());

            int state = 0;
            if ( (*sit)->hasChrgState(state) ) {
                fprintf(out, "N ");
                for (pit = poles.begin(); pit < poles.end(); ++pit) {
                    (*pit)->Charge(state);
                }
                this->CalculateAlphaRigid(poles, out, true);
            }

            state = -1;
            if ( (*sit)->hasChrgState(state) ) {
                fprintf(out, "A ");
                for (pit = poles.begin(); pit < poles.end(); ++pit) {
                    (*pit)->Charge(state);
                }
                this->CalculateAlphaRigid(poles, out, true);
            }

            state = +1;
            if ( (*sit)->hasChrgState(state) ) {
                fprintf(out, "C ");
                for (pit = poles.begin(); pit < poles.end(); ++pit) {
                    (*pit)->Charge(state);
                }
                this->CalculateAlphaRigid(poles, out, true);
            }

            for (pit = poles.begin(); pit < poles.end(); ++pit) {
                (*pit)->Charge(0);
            }

            fprintf(out, " \n");
        }

        return 0;
    }


    // +++++++++++++++++++++++++++++++++ //
    // Create + start threads (Site Ops) //
    // +++++++++++++++++++++++++++++++++ //

    vector<SiteOpQMultipole*> siteOps;

    // Evaluate energies for range of sites?
    if (_listSegIds.size() < 1) {
        _nextSite = top->Segments().begin();

        // Forward to first segment specified in options
        for ( ; (*_nextSite)->getId() != this->_firstSeg &&
                  _nextSite < top->Segments().end(); ++_nextSite) { ; }
    }
    // Evaluate energies for list of sites?
    else {
        for (unsigned int i = 0; i < _listSegIds.size(); ++i) {
            int segId = _listSegIds[i];
            _listSegs.push_back(top->getSegment(segId));
        }
        _nextSite = _listSegs.begin();
    }


    for (unsigned int id = 0; id < _nThreads; id++) {
        SiteOpQMultipole *newOp = new SiteOpQMultipole(id, top, this);
        siteOps.push_back(newOp);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        siteOps[id]->InitSlotData(top);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        siteOps[id]->Start();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        siteOps[id]->WaitDone();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        delete siteOps[id];
    }

    siteOps.clear();


    // +++++++++++++++++++++ //
    // Print numbers to file //
    // +++++++++++++++++++++ //

    if (_energies2File) {

        FILE *out;
        string topOutFile = "Frame"
                          + boost::lexical_cast<string>(top->getDatabaseId())
                          + "_" + _outFile;
        out = fopen(topOutFile.c_str(), "w");

        vector< Segment* > ::iterator sit;
        for (sit = _log_evaluated_segs.begin();
             sit < _log_evaluated_segs.end();
             ++sit) {

            fprintf(out, "%4d ", (*sit)->getId() );
            fprintf(out, "%4s ", (*sit)->getName().c_str() );

            // Energies
            if ((*sit)->hasChrgState(0)) {
                fprintf(out, "   0 %3.8f   ", (*sit)->getEMpoles(0) );
            }
            if ((*sit)->hasChrgState(-1)) {
                fprintf(out, "  -1 %3.8f   ", (*sit)->getEMpoles(-1) );
            }

            if ((*sit)->hasChrgState(+1)) {
                fprintf(out, "  +1 %3.8f   ", (*sit)->getEMpoles(+1) );
            }

            // Iterations
            if ((*sit)->hasChrgState(0)) {
                fprintf(out, "   0 %3d   ", _log_seg_iter[(*sit)->getId()][1]);
            }
            if ((*sit)->hasChrgState(-1)) {
                fprintf(out, "  -1 %3d   ", _log_seg_iter[(*sit)->getId()][0]);
            }

            if ((*sit)->hasChrgState(+1)) {
                fprintf(out, "  +1 %3d   ", _log_seg_iter[(*sit)->getId()][2]);
            }

            // Polarizable sphere
            fprintf(out, "   SPH %4d   ",
                    _log_seg_sphereSize[(*sit)->getId()]);

            fprintf(out, "   %4.7f %4.7f %4.7f   ",
                    _log_seg_com[(*sit)->getId()].getX(),
                    _log_seg_com[(*sit)->getId()].getY(),
                    _log_seg_com[(*sit)->getId()].getZ());

            fprintf(out, " \n");

        }
    }

    this->PostProcess(top);
    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //


/**
 * Provides next segment for thread upon request
 * @return Pointer to next segment
 */
Segment *QMultipole::RequestNextSite(int opId, Topology *top) {

    _nextSiteMutex.Lock();

    Segment *workOnThis;


    // Evaluate energies for range of sites?
    if (_listSegs.size() < 1) {
        if (_nextSite == top->Segments().end()) {
            workOnThis = NULL;
        }
        else if ((*_nextSite)->getId() > _lastSeg && _lastSeg > 0) {
            workOnThis = NULL;
        }
        else {
            workOnThis = *_nextSite;
            _nextSite++;
            cout << endl << "... ... Thread " << opId
                 << " evaluating energy for site "
                 << workOnThis->getId() << flush;
        }
    }

    // Evaluate energies for list of sites?
    else {
        if (_nextSite == _listSegs.end()) {
            workOnThis = NULL;
        }
        else {
            workOnThis = *_nextSite;
            _nextSite++;
            cout << endl << "... ... Thread " << opId
                 << " evaluating energy for site "
                 << workOnThis->getId() << flush;
        }
    }




    _nextSiteMutex.Unlock();

    return workOnThis;
}


// +++++++++++++++++++++++++++++ //
// SiteOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void QMultipole::SiteOpQMultipole::Run(void) {

    while (true) {

        _seg = _master->RequestNextSite(_id, _top);

        if (_seg == NULL) { break; }
        else { this->EvalSite(_top, _seg); }
    }
}


/**
 * Creates private copies of polar sites in topology and
 * arranges them by segments so as to avoid mutexes.
 */
void QMultipole::SiteOpQMultipole::InitSlotData(Topology *top) {


    vector< Segment* > ::iterator sitRef;
    vector< vector<APolarSite*> > ::iterator sitNew;
    vector< APolarSite* > ::iterator pitRef;
    vector< APolarSite* > ::iterator pitNew;

    _APolarSites.resize(top->Segments().size());
    assert(top->Segments().size() == _APolarSites.size());

    for (sitRef = top->Segments().begin(), sitNew = _APolarSites.begin();
         sitRef < top->Segments().end();
         ++sitRef, ++sitNew) {

        (*sitNew).resize((*sitRef)->APolarSites().size());

        for (pitRef = (*sitRef)->APolarSites().begin(),
             pitNew = (*sitNew).begin();
             pitRef < (*sitRef)->APolarSites().end();
             ++pitRef, ++ pitNew) {

            *pitNew = new APolarSite();
            (*pitNew)->ImportFrom(*pitRef, "full");
            (*pitNew)->Charge(0);
        }
    }
}


/**
 * Deletes all private thread data, in particular copies of polar sites
 */
QMultipole::SiteOpQMultipole::~SiteOpQMultipole() {

    vector< vector<APolarSite*> > ::iterator sit;
    vector<APolarSite*> ::iterator pol;

    for (sit = _APolarSites.begin(); sit < _APolarSites.end(); ++sit) {
        for (pol = (*sit).begin(); pol < (*sit).end(); ++pol) {
            delete *pol;
        }
        (*sit).clear();
    }
    _APolarSites.clear();
}


/**
 * Calculate electrostatic + polarization energy for segment
 * ... Determine polarization sphere via cut-off
 * ... Check which charges states are available
 * ... For each charge state: Charge, Induce, Energy, Depolarize
 */
void QMultipole::SiteOpQMultipole::EvalSite(Topology *top, Segment *seg) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    // ++++++++++++++++++++++++++ //
    // Define polarization sphere //
    // ++++++++++++++++++++++++++ //

    this->_segsPolSphere.clear(); // <- Segments    within cutoff
    this->_segsOutSphere.clear(); // <- Segments    within cutoff1, cutoff2
    this->_polsPolSphere.clear(); // <- Polar sites within cutoff
    this->_polsOutSphere.clear(); // <- Polar sites within cutoff1, cutoff2

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

        double r12 = abs(_top->PbShortestConnect((*sit)->getPos(),
                                                    seg->getPos()));

        if      ( r12 > _master->_cutoff2) { continue; }

        else if ( r12 > _master->_cutoff ) {
            _segsOutSphere.push_back(*sit);
            _polsOutSphere.push_back( _APolarSites[(*sit)->getId() - 1] );
        }
        else {
            _segsPolSphere.push_back(*sit);
            _polsPolSphere.push_back( _APolarSites[(*sit)->getId() - 1] );
        }
    }


//    FILE *out;
//    string shellFile = "OuterShell.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    for (sit = _segsOutSphere.begin(); sit < _segsOutSphere.end(); ++sit) {
//        (*sit)->WritePDB(out, "Multipoles", "");
//    }
//    fclose(out);



    // +++++++++++++++++++++++++ //
    // Investigate charge states //
    // +++++++++++++++++++++++++ //

    if (_master->_maverick) {
        cout << endl
             << "... ... ... Segments in polarization sphere: "
             << _segsPolSphere.size()
             << "; segments in static shell: "
             << _segsOutSphere.size()
             << flush;
    }

    this->Depolarize();

    // Make sure, outer shell also depolarized, according to strategy
    vector< vector<APolarSite*> > ::iterator sit1;
    vector<APolarSite*> ::iterator pit;
    for (sit1 = _polsOutSphere.begin(); sit1 < _polsOutSphere.end(); ++sit1) {
        for (pit = (*sit1).begin(); pit < (*sit1).end(); ++pit) {
            (*pit)->Depolarize();
        }
    }

    // Keep track of number of iterations needed until converged
    vector<int> iters;
    iters.resize(3);

    int state = 0;
    if (_seg->hasChrgState(state)) {

        if (_master->_maverick) {
            cout << endl << "... ... ... Seg " << seg->getId()
                         << " Charge " << state << flush;
        }

        this->Charge(state);
        int iter = 0;
        if (_master->_induce) iter = this->Induce(state);
        double EState = 0.0;
        if (_master->_induce) EState = this->Energy(state);
        else                 EState = this->EnergyStatic(state);
        _seg->setEMpoles(state, int2eV * EState);
        iters[state+1] = iter;



//        FILE *out;
//        string toFile = boost::lexical_cast<string>(_seg->getId())
//                      + "_R_absU_N.dat";
//        out = fopen(toFile.c_str(), "w");
//        vector< APolarSite* > ::iterator pit;
//        vector< vector< APolarSite* > > ::iterator sit;
//        vec CoM = _seg->getPos();
//        for (sit = _polsPolSphere.begin(); sit < _polsPolSphere.end(); ++sit) {
//            for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
//                vec r  = _top->PbShortestConnect((*pit)->getPos(), CoM);
//                int segId = (*pit)->getSegment()->getId();
//                double R = abs(r);
//                vec u = vec( (*pit)->U1x, (*pit)->U1y, (*pit)->U1z );
//                // double U = abs(u);
//                fprintf(out, "%4.7f %4.7f %4.7f %4.7f  %4.7f %4.7f %4.7f %4d\n",
//                              R, r.getX(), r.getY(), r.getZ(), u.getX(), u.getY(), u.getZ(), segId);
//            }
//        }
//        fclose(out);


//        // For visual check of convergence
//        string mpNAME = "seg1_cation.dat";
//        FILE *mpDat = NULL;
//        mpDat = fopen(mpNAME.c_str(), "w");
//        vector<APolarSite*>::iterator pit;
//        for (pit = _polsPolSphere[0].begin();
//             pit < _polsPolSphere[0].end();
//             ++pit) {
//            (*pit)->PrintInfoVisual(mpDat);
//        }
//        for (pit = _polsPolSphere[1].begin();
//             pit < _polsPolSphere[1].end();
//             ++pit) {
//            (*pit)->PrintInfoVisual(mpDat);
//        }
//        fclose(mpDat);

        // this->Depolarize(); // OVERRIDE
        vector< vector<APolarSite*> > ::iterator sit1;
        vector< APolarSite* > ::iterator pit1;
        for (sit1 = _polsPolSphere.begin(); sit1 < _polsPolSphere.end(); ++sit1) {
            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                (*pit1)->ResetFieldU();
                (*pit1)->ResetFieldP();
                (*pit1)->ResetU1Hist();
            }
        }
    }

    state = -1;
    if (_seg->hasChrgState(state)) {

        if (_master->_maverick) {
            cout << endl << "... ... ... Seg " << seg->getId()
                         << " Charge " << state << flush;
        }

        this->Charge(state);
        int iter = 0;
        if (_master->_induce) iter = this->Induce(state);
        double EState = 0.0;
        if (_master->_induce) EState = this->Energy(state);
        else                 EState = this->EnergyStatic(state);
        _seg->setEMpoles(state, int2eV * EState);
        iters[state+1] = iter;

        /*
        // For visual check of convergence
        string mpNAME = "seg1_cation.dat";
        FILE *mpDat = NULL;
        mpDat = fopen(mpNAME.c_str(), "w");
        vector<APolarSite*>::iterator pit;
        for (pit = _polsPolSphere[0].begin();
             pit < _polsPolSphere[0].end();
             ++pit) {
            (*pit)->PrintInfoVisual(mpDat);
        }
        fclose(mpDat);
        */

        // this->Depolarize(); // OVERRIDE
        vector< vector<APolarSite*> > ::iterator sit1;
        vector< APolarSite* > ::iterator pit1;
        for (sit1 = _polsPolSphere.begin(); sit1 < _polsPolSphere.end(); ++sit1) {
            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                (*pit1)->ResetFieldU();
                (*pit1)->ResetFieldP();
                (*pit1)->ResetU1Hist();
            }
        }
    }

    state = +1;
    if (_seg->hasChrgState(state)) {

        if (_master->_maverick) {
            cout << endl << "... ... ... Seg " << seg->getId()
                         << " Charge " << state << flush;
        }

        this->Charge(state);
        int iter = 0;
        if (_master->_induce) iter = this->Induce(state);
        double EState = 0.0;
        if (_master->_induce) EState = this->Energy(state);
        else                  EState = this->EnergyStatic(state);
        _seg->setEMpoles(state, int2eV * EState);
        iters[state+1] = iter;



//        FILE *out;
//        string toFile = boost::lexical_cast<string>(_seg->getId())
//                      + "_R_absU_C.dat";
//        out = fopen(toFile.c_str(), "w");
//        vector< APolarSite* > ::iterator pit;
//        vector< vector< APolarSite* > > ::iterator sit;
//        vec CoM = _seg->getPos();
//        for (sit = _polsPolSphere.begin(); sit < _polsPolSphere.end(); ++sit) {
//            for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
//                vec r  = _top->PbShortestConnect((*pit)->getPos(), CoM);
//                int segId = (*pit)->getSegment()->getId();
//                double R = abs(r);
//                vec u = vec( (*pit)->U1x, (*pit)->U1y, (*pit)->U1z );
//                // double U = abs(u);
//                fprintf(out, "%4.7f %4.7f %4.7f %4.7f  %4.7f %4.7f %4.7f %4d \n",
//                              R, r.getX(), r.getY(), r.getZ(), u.getX(), u.getY(), u.getZ(), segId);
//            }
//        }
//        fclose(out);


//        // For visual check of convergence
//        string mpNAME = "seg1_neutral.dat";
//        FILE *mpDat = NULL;
//        mpDat = fopen(mpNAME.c_str(), "w");
//        vector<APolarSite*>::iterator pit;
//        for (pit = _polsPolSphere[0].begin();
//             pit < _polsPolSphere[0].end();
//             ++pit) {
//            (*pit)->PrintInfoVisual(mpDat);
//        }
//        for (pit = _polsPolSphere[1].begin();
//             pit < _polsPolSphere[1].end();
//             ++pit) {
//            (*pit)->PrintInfoVisual(mpDat);
//        }
//        fclose(mpDat);


        this->Depolarize();

    }

    this->_master->_logMutex.Lock();
    _master->_log_seg_iter[_seg->getId()] = iters;
    _master->_log_seg_sphereSize[_seg->getId()] = _polsPolSphere.size();
    _master->_log_seg_com[_seg->getId()] = _seg->getPos();
    _master->_log_evaluated_segs.push_back(_seg);
    this->_master->_logMutex.Unlock();

    this->Charge(0);

}


/**
 * Charges polar sites in segment to state specified
 * @param state
 */
void QMultipole::SiteOpQMultipole::Charge(int state) {

    vector< APolarSite* > ::iterator pit;
    for (pit = _APolarSites[_seg->getId()-1].begin();
         pit < _APolarSites[_seg->getId()-1].end();
         ++pit) {

        (*pit)->Charge(state);
    }
}


/**
 * Calculates induced dipole moments using the Thole Model + SOR
 */
int QMultipole::SiteOpQMultipole::Induce(int state) {

    double wSOR = (state == 0) ? _master->_wSOR_N : _master->_wSOR_C;
    double eTOL = this->_master->_epsTol;
    int    maxI = this->_master->_maxIter;

    vector< vector<APolarSite*> > ::iterator sit1;
    vector< vector<APolarSite*> > ::iterator sit2;
    vector< APolarSite* > ::iterator pit1;
    vector< APolarSite* > ::iterator pit2;

    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site fields (arising from perm. m'poles) //
    // ++++++++++++++++++++++++++++++++++++++++++++++ //

//    cout << endl << "... ... ... 0th-order field" << flush;
    for (sit1 = _polsPolSphere.begin();
         sit1 < _polsPolSphere.end();
         ++sit1) {
    for (sit2 = sit1 + 1;
         sit2 < _polsPolSphere.end();
         ++sit2) {

         for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
         for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

             _actor.FieldPerm(*(*pit1), *(*pit2));
         }}
    }}

    // +++++++++++++++++++ //
    // 1st-order induction //
    // +++++++++++++++++++ //

//    cout << " | Induce " << endl;
    if (state == 0) { // OVERRIDE
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

             for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                 (*pit1)->InduceDirect();
             }
        }
    }



    // ++++++++++++++++++++++ //
    // Higher-order induction //
    // ++++++++++++++++++++++ //

    int iter = 0;
//    boost::progress_timer T;
    for ( ; iter < maxI; ++iter) {

        // Reset fields FUx, FUy, FUz
//        cout << "\r... ... ... Reset (" << iter << ")" << flush;
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                (*pit1)->ResetFieldU();
            }
        }

        // Intra-site contribution to induction field
//        cout << " | Intra-Site (" << iter << ")" << flush;
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
            for (pit2 = pit1 + 1;        pit2 < (*sit1).end(); ++pit2) {

                _actor.FieldIndu(*(*pit1),*(*pit2));
            }}
        }

        // Inter-site contribution to induction field
//        cout << " | Inter-Site (" << iter << ")" << flush;
        //boost::progress_display show_progress( _polsPolSphere.size() );
//        T.restart();
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {
        for (sit2 = sit1 + 1;
             sit2 < _polsPolSphere.end();
             ++sit2) {

            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
            for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

                _actor.FieldIndu(*(*pit1), *(*pit2));
            }}
        }}
//        cout << " | dt " << T.elapsed() << flush;

        // Induce again
//        cout << " | Induce (" << iter << ")" << flush;
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

             for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                 (*pit1)->Induce(wSOR);
             }
        }

        // Check for convergence
//        cout << " | Check (" << iter << ")" << flush;
        bool converged = true;
        double maxdU = -1;
        double avgdU = 0.0;
        int    baseN = 0;
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

             for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                 double dU = (*pit1)->HistdU();
                 avgdU += dU;
                 ++baseN;
                 if ( dU > maxdU ) { maxdU = dU; }
                 if ( dU > eTOL ) { converged = false; }
             }
        }
        avgdU /= baseN;
        if (avgdU < eTOL/10.) { converged = true; }

//        cout << " | MAX dU " << maxdU
//             << " | AVG dU " << avgdU
//             << " | SOR " << wSOR << flush;

        // Break if converged
        if      (converged) {
            break;
        }
        else if (iter == maxI - 1) {
            this->_master->LockCout();
            cout << endl << "... ... ... WARNING Induced multipoles for site "
                 << _seg->getId() << " - state " << state
                 << " did not converge to precision: "
                 << " AVG dU:U " << avgdU << flush;
            this->_master->UnlockCout();
            break;
        }
    }

    return iter;

//    cout << endl << "... ... ... State " << state
//         << " - wSOR " << wSOR
//         << " - Iterations " << iter << flush;
}


double QMultipole::SiteOpQMultipole::EnergyStatic(int state) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    _actor.ResetEnergy();
    double E_Tot = 0.0;

    vector< Segment* > ::iterator seg1;
    vector< Segment* > ::iterator seg2;
    vector< vector<APolarSite*> > ::iterator sit1;
    vector< vector<APolarSite*> > ::iterator sit2;
    vector< APolarSite* > ::iterator pit1;
    vector< APolarSite* > ::iterator pit2;

    vector< APolarSite* > central = _APolarSites[ _seg->getId() - 1 ];


    for (seg1 = _segsPolSphere.begin(); seg1 < _segsPolSphere.end(); ++seg1) {

        int id = (*seg1)->getId();

        if (id == _seg->getId()) {
            continue;
        }

        for (pit1 = _APolarSites[id-1].begin();
             pit1 < _APolarSites[id-1].end();
             ++pit1) {
        for (pit2 = central.begin();
             pit2 < central.end();
             ++pit2) {

             E_Tot += _actor.EnergyInter(*(*pit1), *(*pit2));
        }}
    }

    if (_master->_maverick) {
        cout << endl << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
             << " = (P ~) " << _actor.getEP() * int2eV
             << " + (U ~) " << _actor.getEU_INTER() * int2eV
             << " + (U o) " << _actor.getEU_INTRA() * int2eV
             << " , statics only. "
             << flush;
    }

    return E_Tot;
}


/**
 * Calculates electrostatic + induction energy of segment up to Q2-Q2
 */
double QMultipole::SiteOpQMultipole::Energy(int state) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    _actor.ResetEnergy();
    double E_Tot = 0.0;

    vector< Segment* > ::iterator seg1;
    vector< Segment* > ::iterator seg2;
    vector< vector<APolarSite*> > ::iterator sit1;
    vector< vector<APolarSite*> > ::iterator sit2;
    vector< APolarSite* > ::iterator pit1;
    vector< APolarSite* > ::iterator pit2;

    // +++++++++++++++++ //
    // Inter-site energy //
    // +++++++++++++++++ //

    for (sit1 = _polsPolSphere.begin(), seg1 = _segsPolSphere.begin();
         sit1 < _polsPolSphere.end();
         ++sit1, ++seg1) {
    for (sit2 = sit1 + 1, seg2 = seg1 + 1;
         sit2 < _polsPolSphere.end();
         ++sit2, ++seg2) {

        for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
        for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

            E_Tot += _actor.E_f(*(*pit1), *(*pit2));
            E_Tot += _actor.E_m(*(*pit1), *(*pit2));
            E_Tot += _actor.E_m(*(*pit2), *(*pit1));

        }}
    }}

    // +++++++++++++++++ //
    // Intra-site energy //
    // +++++++++++++++++ //

//    // Not needed here
//    for (sit1 = _polsPolSphere.begin(), seg1 = _segsPolSphere.begin();
//         sit1 < _polsPolSphere.end();
//         ++sit1, ++seg1) {
//
//        for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
//        for (pit2 = pit1 + 1; pit2 < (*sit1).end(); ++pit2) {
//
//            E_Tot += _actor.EnergyIntra(*(*pit1), *(*pit2)); // OVERRIDE
//
//        }}
//    }

    // ++++++++++++++++++ //
    // Outer-Shell energy //
    // ++++++++++++++++++ //

    double E_Out = 0.0;

    vector< APolarSite* > central = _APolarSites[ _seg->getId() - 1 ];

    for (sit1 = _polsOutSphere.begin(); sit1 < _polsOutSphere.end(); ++sit1) {
        for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
            for (pit2 = central.begin(); pit2 < central.end(); ++pit2) {
                E_Out += _actor.E_f(*(*pit1), *(*pit2));
                E_Out += _actor.E_m(*(*pit2), *(*pit1));
            }
        }
    }
    E_Tot += E_Out;


    if (_master->_maverick) {
        cout << endl << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
             << " = (P ~) " << _actor.getEP()       * int2eV
             << " + (U ~) " << _actor.getEU_INTER() * int2eV
             << " + (U o) " << _actor.getEU_INTRA() * int2eV
             << " , (O ~) " << E_Out * int2eV
             << flush;
    }
    
    if (_master->_maverick) {
        cout << endl << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
             << " = (PP ~) " << _actor.getEPP() * int2eV
             << " + (PU ~) " << _actor.getEPU() * int2eV
             << " + (UU ~) " << _actor.getEUU() * int2eV
             << flush;
    }

    return E_Tot;
}


/**
 * Zeroes out induced moments and induction fields for all polar sites
 * in polarizable sphere
 */
void QMultipole::SiteOpQMultipole::Depolarize() {

    vector< vector<APolarSite*> > ::iterator sit;
    vector< APolarSite* > ::iterator pit;

    for (sit = _polsPolSphere.begin(); sit < _polsPolSphere.end(); ++sit) {
        for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {

            (*pit)->Depolarize();
        }
    }
}


// +++++++++++++++++++++++++++ //
// Interactor Member Functions //
// +++++++++++++++++++++++++++ //

/**
 * Used in ESP calculator (initialize stage of EMultipole)
 */
inline double QMultipole::Interactor::PotentialPerm(vec r,
                                                     APolarSite &pol) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol.getPos() - r;
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;




    rbx = - pol._locX * e12;
    rby = - pol._locY * e12;
    rbz = - pol._locZ * e12;

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
 * Used in ESF calculator (initialize stage of EMultipole)
 */
inline vec QMultipole::Interactor::FieldPermESF(vec r,
                                                 APolarSite &pol) {
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
inline void QMultipole::Interactor::FieldInduAlpha(APolarSite &pol1,
                                                    APolarSite &pol2) {
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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
inline void QMultipole::Interactor::FieldIndu(APolarSite &pol1,
                                               APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //          This implies that induced = - alpha * field
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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
inline void QMultipole::Interactor::FieldPerm(APolarSite &pol1,
                                               APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    //          This implies that induced = - alpha * field
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
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
inline double QMultipole::Interactor::EnergyIntra(APolarSite &pol1,
                                                   APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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
inline double QMultipole::Interactor::EnergyInter(APolarSite &pol1,
                                                   APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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
inline double QMultipole::Interactor::EnergyInterESP(APolarSite &pol1,
                                                      APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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


inline double QMultipole::Interactor::E_f(APolarSite &pol1,
                                          APolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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
        
    EP += epp;
    EU_INTER += 0.5 * epu;
        
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
    

    EPP += epp;
    EPU += epu;
    EUU += euu;
    return epp + epu + euu;
}

inline double QMultipole::Interactor::E_m(APolarSite &pol1,
                                          APolarSite &pol2) {
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _em->_aDamp;
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

    double epp = 0.0; // <- Interaction perm. <> perm.
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
    epp *= -0.5;
    epu *= -0.5;

    EPP += epp;
    EPU += epu;
    EUU += euu;
    return epp + epu + euu;
}



}}







#endif
