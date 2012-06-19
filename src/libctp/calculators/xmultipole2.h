#ifndef XMULTIPOLE2_H
#define XMULTIPOLE2_H


#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

class XMP : public QMCalculator
{

public:

    XMP() {};
   ~XMP() {};

    string              Identify() { return "XMultipole"; }
    void                Initialize(Topology *, Property *);
    bool                EvaluateFrame(Topology *top);

    // +++++++++++++++++++++++++++++++ //
    // MULTIPOLE ALLOCATION FROM INPUT //
    // +++++++++++++++++++++++++++++++ //

    void                Collect_JOB(string job_file, Topology *top);
    void                Collect_EMP(string emp_file, Topology *top);
    void                Collect_MPS(Topology *top);
    void                Collect_XML(string xml_file, Topology *top);

    void                Create_MPOLS(Topology *top);
    vector<PolarSite*>  Map_MPols_To_Seg(vector<PolarSite*> &, Segment *);
    vector<PolarSite*>  Parse_GDMA(string mps_file, int state);

    // +++++++++++++++++++++++++++ //
    // Multipole Interaction Class //
    // +++++++++++++++++++++++++++ //

    class XInteractor
    {
    public:

        XInteractor(Topology *top, XMP *xm) : _top(top), _xm(xm) {};
        XInteractor() : _top(NULL), _xm(NULL) {};
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


        inline double EnergyInter(PolarSite &pol1, PolarSite &pol2);
        inline double EnergyInterESP(PolarSite &pol1, PolarSite &pol2);
        inline double EnergyIntra(PolarSite &pol1, PolarSite &pol2);
        inline void FieldPerm(PolarSite &pol1, PolarSite &pol2);
        inline vec  FieldPermESF(vec r, PolarSite &pol);
        inline void FieldIndu(PolarSite &pol1, PolarSite &pol2);
        inline void FieldInduAlpha(PolarSite &pol1, PolarSite &pol2);

        inline double PotentialPerm(vec r, PolarSite &pol);

        void ResetEnergy() { EP = EU_INTER = EU_INTRA = 0.0; }
        double &getEP() { return EP; }
        double &getEU_INTER() { return EU_INTER; }
        double &getEU_INTRA() { return EU_INTRA; }

    private:

        double EP;       //   <- Interaction permanent multipoles (inter-site)
        double EU_INTRA; //   <- Interaction induction multipoles (intra-site)
        double EU_INTER; //   <- Interaction induction multipoles (inter-site)


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
        XMP             *_xm;

    };

    // +++++++++++++++++++++++++++ //
    // Job Operator (Thread class) //
    // +++++++++++++++++++++++++++ //

    class JobXMP : public Thread
    {
    public:

        JobXMP(int id,     Topology *top,        XMP *master)
                : _id(id),      _top(top),    _master(master)
                   { _actor = XInteractor(top,_master); };

       ~JobXMP();

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void InitSlotData(Topology *top) { ; }
        void Run(void) { ; }
        void EvalSite(Topology *top, Segment *seg) { ; }



    public:

        int                           _id;
        Topology                     *_top;
        XMP                          *_master;

        vector< Segment* >           _segsPolSphere; // Segments    in c/o 0-1
        vector< Segment* >           _segsOutSphere; // Segments    in c/0 1-2
        vector< vector<PolarSite*> > _polsPolSphere; // Polar sites in c/o 0-1
        vector< vector<PolarSite*> > _polsOutSphere; // Polar sites in c/o 1-2
        vector< vector<PolarSite*> > _polarSites;    // Copy of top polar sites
        XInteractor                  _actor;
    };


private:

    // ++++++++++++++++++++++++++++++++++++++++ //
    // MULTIPOLE ALLOCATION TO SEGMENTS / PAIRS //
    // ++++++++++++++++++++++++++++++++++++++++ //

    string                         _job_file;
    string                         _emp_file;
    string                         _xml_file;

    // Job info : JOB_ID PAIR_ID MPS_1 MPS_2 TAG
    map<int,string>                _jobId_jobTag;
    map<int,int>                   _jobId_pairId;
    map<int,string>                _jobId_mpsFile1;
    map<int,string>                _jobId_mpsFile2;
    map<int, pair<int,int> >       _pairId_seg1Id_seg2Id;

    // Allocate .mps files to segments (neutral, electron, hole)
    map<int,string>                 _segId_mpsFile_n;   // <=> charge state  0
    map<int,string>                 _segId_mpsFile_e;   // <=> charge state -1
    map<int,string>                 _segId_mpsFile_h;   // <=> charge state +1  

    // Store polar site containers, one for each GDMA (= .mps) file
    map<string,vector<PolarSite*> > _mpsFile_pSites;

    // Allocate polar sites to fragments in segments
    map<string, bool>               _map2md;
    map<string, vector<int> >       _alloc_frag_mpoleIdx;
    map<string, vector<string> >    _alloc_frag_mpoleName;
    map<string, vector<int> >       _alloc_frag_trihedron;
    map<string, vector<double> >    _alloc_frag_weights;

    string                          _pdb_check;

};

void XMP::Initialize(Topology *top, Property *options) {

    string key  = "options.xmultipole.control";
    _job_file   = options->get(key+".job_file").as<string>();
    _emp_file   = options->get(key+".emp_file").as<string>();

    if (options->exists(key+".pdb_check")) {
        _pdb_check = options->get(key+".pdb_check").as<string>();
    }
    else { _pdb_check = ""; }

    key         = "options.xmultipole.multipoles";
    _xml_file   = options->get(key).as<string>();

}


void XMP::Collect_XML(string xml_file, Topology *top) {

    cout << endl << "... ... Allocate polar sites to fragments. " << flush;

    string allocFile = xml_file;

    // ++++++++++++++++++++++++++++++++ //
    // Load polar-site indices from XML //
    // ++++++++++++++++++++++++++++++++ //

    // => Output to maps:
    map<string, vector<int> > alloc_frag_mpoleIdx;
    map<string, vector<string> > alloc_frag_mpoleName;
    map<string, vector<int> > alloc_frag_trihedron;
    map<string, vector<double> > alloc_frag_weights;

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


    string key = "topology.molecules.molecule";
    list<Property *> mols = allocation.Select(key);
    list<Property *>::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); molit++) {

        string molName = (*molit)->get("name").as<string> ();

        key = "segments.segment";
        list<Property *> segs = (*molit)->Select(key);
        list<Property *>::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); segit++) {

            string segName = (*segit)->get("name").as<string> ();

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

                vector<int> trihedron = (*fragit)->get("localframe")
                                                .as< vector<int> >();

                vector<double> weights = (*fragit)->get("weights")
                                                .as< vector<double> >();

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
                alloc_frag_trihedron[mapKeyName] = trihedron;
                alloc_frag_weights[mapKeyName] = weights;
            }
        }
    }

    _alloc_frag_mpoleIdx    = alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName   = alloc_frag_mpoleName;
    _alloc_frag_trihedron   = alloc_frag_trihedron;
    _alloc_frag_weights     = alloc_frag_weights;
}


void XMP::Collect_JOB(string job_file, Topology *top) {

    QMNBList &nblist = top->NBList();

    std::string line;
    std::ifstream intt;
    intt.open(job_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

// Sample line
// # JOB_ID TAG  PAIR_ID SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS
//   1      E_CT 3819    182     C60       c60.mps  392     DCV       dcv.mps

            int jobId       = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            int pairId      = boost::lexical_cast<int>(split[2]);

            int seg1Id      = boost::lexical_cast<int>(split[3]);
            string seg1Name = split[4];
            string seg1mps  = split[5];

            int seg2Id      = boost::lexical_cast<int>(split[6]);
            string seg2Name = split[7];
            string seg2mps  = split[8];

            Segment *seg1   = top->getSegment(seg1Id);
            Segment *seg2   = top->getSegment(seg2Id);
            QMPair  *qmp    = nblist.FindPair(seg1,seg2);

            if (qmp == NULL) {
                cout << endl << "ERROR: '" << job_file
                     << "': No such pair "
                     << pairId << " " << seg1Name << " " << seg2Name << ". "
                     << flush;
                throw runtime_error("Pair specs do not match topology.");
            }

            if ( _jobId_jobTag.count(jobId) ) {
                cout << endl << "ERROR: '" << job_file
                     << "': Job-ID " << jobId << " exists more "
                        "than once. Abort." << endl;
                throw runtime_error("Rework job file.");
            }

            _jobId_jobTag[jobId]            = tag;
            _jobId_pairId[jobId]            = pairId;
            _jobId_mpsFile1[jobId]          = seg1mps;
            _jobId_mpsFile2[jobId]          = seg2mps;
            _pairId_seg1Id_seg2Id[pairId]   = pair<int,int>(seg1Id,seg2Id);


        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl << "... ... Registered " << _jobId_jobTag.size() << " jobs. "
         << flush;

}


void XMP::Collect_MPS(Topology *top) {

    // +++++++++++++ //
    // Parse + Store //
    // +++++++++++++ //

    map<int,string> ::iterator misit;

    // Neutral species
    for (misit = _segId_mpsFile_n.begin();
         misit != _segId_mpsFile_n.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, 0); }
    }

    // Anion
    for (misit = _segId_mpsFile_e.begin();
         misit != _segId_mpsFile_e.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, -1); }
    }

    // Cation
    for (misit = _segId_mpsFile_h.begin();
         misit != _segId_mpsFile_h.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    // Job Seg1
    for (misit = _jobId_mpsFile1.begin();
         misit != _jobId_mpsFile1.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    // Job Seg2
    for (misit = _jobId_mpsFile2.begin();
         misit != _jobId_mpsFile2.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    cout << endl << "... ... Parsed a total of " << _mpsFile_pSites.size()
            << " multipole-definition (.mps) files. " << flush;
}


void XMP::Collect_EMP(string emp_file, Topology *top) {

    std::string line;
    std::ifstream intt;
    intt.open(emp_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

            // Line sample:
            // 1022 C60 C60_n.mps C60_e.mps C60_h.mps

            int segId       = boost::lexical_cast<int>(split[0]);
            string segName  = boost::lexical_cast<string>(split[1]);

            // Some compliance checks
            assert( split.size() == 5 );                                        // <- File format correct?
            assert( top->getSegment(segId)->getName() == segName );             // <- Input matches topology?

            _segId_mpsFile_n[segId] = split[2];
            _segId_mpsFile_e[segId] = split[3];
            _segId_mpsFile_h[segId] = split[4];

        } // Exit loop over lines
    }
    else { cout << endl << "ERROR: No such file " << emp_file << endl; }

    assert( _segId_mpsFile_n.size() == top->Segments().size() );                // <- Input for all segments?

}


void XMP::Create_MPOLS(Topology *top) {

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg                = *sit;
        int segId                   = seg->getId();

        bool map2md                 = _map2md[seg->getName()];

        string mps_n                = _segId_mpsFile_n[segId];
        string mps_h                = _segId_mpsFile_h[segId];
        string mps_e                = _segId_mpsFile_e[segId];

        vector<PolarSite*> pols_n   = _mpsFile_pSites[mps_n];
        vector<PolarSite*> pols_h   = _mpsFile_pSites[mps_h];
        vector<PolarSite*> pols_e   = _mpsFile_pSites[mps_e];

        // Merge polar sites
        assert(pols_n.size() == pols_h.size());
        assert(pols_n.size() == pols_e.size());

        for (int i = 0; i < pols_n.size(); i++) {
            pols_n[i]->setQs( pols_h[i]->getQs(+1), +1 );
            pols_n[i]->setPs( pols_h[i]->getPs(+1), +1 );
        }
        for (int i = 0; i < pols_n.size(); ++i) {
            pols_n[i]->setQs(pols_e[i]->getQs(-1), -1 );
            pols_n[i]->setPs(pols_e[i]->getPs(-1), -1 );
        }

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            // Retrieve polar-site data for this fragment
            string idkey                = frag->getName() + seg->getName()
                                        + seg->getMolecule()->getName();
            vector<int> polesInFrag     = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag  = _alloc_frag_mpoleName.at(idkey);
            vector<double> weightsInFrag= _alloc_frag_weights.at(idkey);

            if (map2md && polesInFrag.size() != frag->Atoms().size()) {
                cout << endl
                     << "ERROR: Segment " << seg->getName()
                     << " Fragment " << frag->getName()
                     << ": MAP2MD = TRUE requires same number of polar "
                     << "sites as there are atoms to perform mapping. "
                     << endl;
                throw runtime_error("Check mapping or switch to map2md = 0");
            }

            matrix rotateMP2MD;
            vec translateMP2MD;

            // Determine transformation matrices to execute mapping
            if (!map2md) {

                vector<PolarSite*> trihedron_pol;
                vector<Atom*>      trihedron_atm;

                vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
                vector<int> ::iterator iit;
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    trihedron_pol.push_back(pols_n[(*iit)-1]);

                    // Find fragment-internal position = id of leg
                    for (int i = 0; i < polesInFrag.size(); ++i) {
                        if (polesInFrag[i] == (*iit)) {
                            trihedron_atm.push_back(frag->Atoms()[i]);                            
                        }
                    }
                }


                int symmetry = trihedron_pol.size();
                assert (trihedron_pol.size() == trihedron_atm.size() );                

                vec xMD, yMD, zMD;
                vec xQM, yQM, zQM;

                if (symmetry == 3) {
                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r3MD = trihedron_atm[2]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();
                    vec r3QM = trihedron_pol[2]->getPos();

                    xMD = r2MD - r1MD;
                    yMD = r3MD - r1MD;
                    xQM = r2QM - r1QM;
                    yQM = r3QM - r1QM;

                    zMD = xMD ^ yMD;
                    zQM = xQM ^ yQM;

                    yMD = zMD ^ xMD;
                    yQM = zQM ^ xQM;

                    xMD = xMD.normalize();
                    yMD = yMD.normalize();
                    zMD = zMD.normalize();
                    xQM = xQM.normalize();
                    yQM = yQM.normalize();
                    zQM = zQM.normalize();
                }

                else if (symmetry = 2) {

                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();

                    xMD = r2MD - r1MD;
                    xQM = r2QM - r1QM;

                    // Normalising not necessary, but recommendable, when doing...
                    xMD = xMD.normalize();
                    xQM = xQM.normalize();

                    vec yMDtmp = vec(0,0,0);
                    vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
       // All components of xMD are unequal to zero => division is safe.
       // Choose vector from plane with xMD * inPlaneVec = 0:
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
       yMDtmp = vec(tmp_x, tmp_y, tmp_z);
       yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
       yQMtmp = vec(tmp_x, tmp_y, tmp_z);
       yQMtmp.normalize();
    }

                    // Now proceed as for symmetry 3
                    zMD = xMD ^ yMDtmp;
                    yMD = zMD ^ xMD;
                    zQM = xQM ^ yQMtmp;
                    yQM = zQM ^ xQM;

                    xMD.normalize();
                    yMD.normalize();
                    zMD.normalize();
                    xQM.normalize();
                    yQM.normalize();
                    zQM.normalize();
                }

                else if (symmetry = 1) {

                    cout << endl
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << endl;

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                else {
                    cout << endl
                         << "NOTE: Invalid definition of local frame in fragment "
                         << frag->getName();
                    cout << ". Assuming point particle for mapping. "
                         << endl;
                    cout << endl
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << endl;

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                matrix rotMD = matrix(xMD, yMD, zMD);
                matrix rotMP = matrix(xQM, yQM, zQM);
                
                rotateMP2MD = rotMD * rotMP.Transpose();


                // ++++++++++++++++++ //
                // Transform fragment //
                // ++++++++++++++++++ //

                
                vec CoMP = vec(0.,0.,0.);
                double W = 0.0;
                for (int i = 0; i < polesInFrag.size(); ++i) {

                    double weight = weightsInFrag[i];
                    
                    vec pos = pols_n[polesInFrag[i]-1]->getPos();
                    
                    CoMP += weight*pos;
                    W += weight;

                }
                CoMP /= W;

                translateMP2MD = frag->getCoMD() - CoMP;

            }            

            // Create polar sites 
            for (int i = 0; i < polesInFrag.size(); i++) {

                string name             = namesInFrag[i];
                int poleId              = polesInFrag[i];

                PolarSite *templ        = pols_n[poleId-1];
                PolarSite *newSite      = top->AddPolarSite(name);
                newSite->ImportFrom(templ);
                seg->AddPolarSite(newSite);
                frag->AddPolarSite(newSite);

                // Shift + rotate
                if (!map2md) {
                    newSite->Translate(translateMP2MD);
                    newSite->Rotate(rotateMP2MD, frag->getCoMD());
                }
                else {
                    vec mdpos = frag->Atoms()[i]->getPos();
                    newSite->setPos(mdpos);
                    if (newSite->getRank() > 0) {
                        cout << endl 
                             << "ERROR: MAP2MD = TRUE prevents use of "
                             << "higher-rank multipoles. "
                             << endl;
                        throw runtime_error("User not paying attention. ");
                    }
                }
            }
        }
    }

    top->setIsEStatified(true);

}


vector<PolarSite*> XMP::Map_MPols_To_Seg(vector<PolarSite*> &pols_n, Segment *seg) {

    vector<PolarSite*> return_pols;
    return_pols.resize(pols_n.size());

    int segId                   = seg->getId();
    bool map2md                 = _map2md[seg->getName()];

    vector<Fragment*> ::iterator fit;
    for (fit = seg->Fragments().begin();
         fit < seg->Fragments().end();
         ++fit) {

        Fragment *frag = *fit;

        // Retrieve polar-site data for this fragment
        string idkey                = frag->getName() + seg->getName()
                                    + seg->getMolecule()->getName();
        vector<int> polesInFrag     = _alloc_frag_mpoleIdx.at(idkey);
        vector<string> namesInFrag  = _alloc_frag_mpoleName.at(idkey);
        vector<double> weightsInFrag= _alloc_frag_weights.at(idkey);

        if (map2md && polesInFrag.size() != frag->Atoms().size()) {
            cout << endl
                 << "ERROR: Segment " << seg->getName()
                 << " Fragment " << frag->getName()
                 << ": MAP2MD = TRUE requires same number of polar "
                 << "sites as there are atoms to perform mapping. "
                 << endl;
            throw runtime_error("Check mapping or switch to map2md = 0");
        }

        matrix rotateMP2MD;
        vec translateMP2MD;

        // Determine transformation matrices to execute mapping
        if (!map2md) {

            vector<PolarSite*> trihedron_pol;
            vector<Atom*>      trihedron_atm;

            vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
            vector<int> ::iterator iit;
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                trihedron_pol.push_back(pols_n[(*iit)-1]);

                // Find fragment-internal position = id of leg
                for (int i = 0; i < polesInFrag.size(); ++i) {
                    if (polesInFrag[i] == (*iit)) {
                        trihedron_atm.push_back(frag->Atoms()[i]);
                    }
                }
            }


            int symmetry = trihedron_pol.size();
            assert (trihedron_pol.size() == trihedron_atm.size() );

            vec xMD, yMD, zMD;
            vec xQM, yQM, zQM;

            if (symmetry == 3) {
                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r3MD = trihedron_atm[2]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();
                vec r3QM = trihedron_pol[2]->getPos();

                xMD = r2MD - r1MD;
                yMD = r3MD - r1MD;
                xQM = r2QM - r1QM;
                yQM = r3QM - r1QM;

                zMD = xMD ^ yMD;
                zQM = xQM ^ yQM;

                yMD = zMD ^ xMD;
                yQM = zQM ^ xQM;

                xMD = xMD.normalize();
                yMD = yMD.normalize();
                zMD = zMD.normalize();
                xQM = xQM.normalize();
                yQM = yQM.normalize();
                zQM = zQM.normalize();
            }

            else if (symmetry = 2) {

                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();

                xMD = r2MD - r1MD;
                xQM = r2QM - r1QM;

                // Normalising not necessary, but recommendable, when doing...
                xMD = xMD.normalize();
                xQM = xQM.normalize();

                vec yMDtmp = vec(0,0,0);
                vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
    // All components of xMD are unequal to zero => division is safe.
    // Choose vector from plane with xMD * inPlaneVec = 0:
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
    yMDtmp = vec(tmp_x, tmp_y, tmp_z);
    yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
    yQMtmp = vec(tmp_x, tmp_y, tmp_z);
    yQMtmp.normalize();
    }

                // Now proceed as for symmetry 3
                zMD = xMD ^ yMDtmp;
                yMD = zMD ^ xMD;
                zQM = xQM ^ yQMtmp;
                yQM = zQM ^ xQM;

                xMD.normalize();
                yMD.normalize();
                zMD.normalize();
                xQM.normalize();
                yQM.normalize();
                zQM.normalize();
            }

            else if (symmetry = 1) {

                cout << endl
                     << "WARNING: Symmetry = 1 for fragment "
                     << frag->getName() << ": This will generate artifacts "
                     << "when mapping higher-rank multipoles (dipoles, ..)."
                     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            else {
                cout << endl
                     << "NOTE: Invalid definition of local frame in fragment "
                     << frag->getName();
                cout << ". Assuming point particle for mapping. "
                     << endl;
                cout << endl
                     << "WARNING: Symmetry = 1 for fragment "
                     << frag->getName() << ": This will generate artifacts "
                     << "when mapping higher-rank multipoles (dipoles, ..)."
                     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            matrix rotMD = matrix(xMD, yMD, zMD);
            matrix rotMP = matrix(xQM, yQM, zQM);

            rotateMP2MD = rotMD * rotMP.Transpose();


            // ++++++++++++++++++ //
            // Transform fragment //
            // ++++++++++++++++++ //


            vec CoMP = vec(0.,0.,0.);
            double W = 0.0;
            for (int i = 0; i < polesInFrag.size(); ++i) {

                double weight = weightsInFrag[i];

                vec pos = pols_n[polesInFrag[i]-1]->getPos();

                CoMP += weight*pos;
                W += weight;

            }
            CoMP /= W;

            translateMP2MD = frag->getCoMD() - CoMP;

        }

        // Create polar sites
        for (int i = 0; i < polesInFrag.size(); i++) {

            string name             = namesInFrag[i];
            int poleId              = polesInFrag[i];

            PolarSite *templ        = pols_n[poleId-1];
            PolarSite *newSite      = new PolarSite(-1, name);
            newSite->ImportFrom(templ);

            // Shift + rotate
            if (!map2md) {
                newSite->Translate(translateMP2MD);
                newSite->Rotate(rotateMP2MD, frag->getCoMD());
            }
            else {
                vec mdpos = frag->Atoms()[i]->getPos();
                newSite->setPos(mdpos);
                if (newSite->getRank() > 0) {
                    cout << endl
                         << "ERROR: MAP2MD = TRUE prevents use of "
                         << "higher-rank multipoles. "
                         << endl;
                    throw runtime_error("User not paying attention. ");
                }
            }

            return_pols.push_back(newSite);

        }
    } // End loop over fragments

    return return_pols;
}


vector<PolarSite*> XMP::Parse_GDMA(string filename, int state) {

    int poleCount = 1;
    double Q0_total = 0.0;
    string units = "";
    bool useDefaultPs = true;

    vector<PolarSite*> poles;
    PolarSite *thisPole = NULL;

    vector<double> Qs; // <- multipole moments
    double         P1; // <- dipole polarizability

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

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
    // ! DCV2T opt
    // ! SP        RB3LYP          6-311+G(d,p)
    // Units bohr
    //
    // C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
    //  -0.3853409355
    //  -0.0002321905   0.2401559510   0.6602334308
    //  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813  -0.50937399
    //  P 1.75


            // Units used
            if ( split[0] == "Units") {
                units = split[1];
                if (units != "bohr" && units != "angstrom") {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }
            }

            // element,  position,  rank limit
            else if ( split.size() == 6 ) {

                Qs.clear();
                P1 = -1.;

                int id = poleCount++;  // <- starts from 1
                string name = split[0];

                double BOHR2NM = 0.0529189379;
                double ANGSTROM2NM = 0.1;
                double x, y, z;

                if (units == "bohr") {
                    x = BOHR2NM * boost::lexical_cast<double>(split[1]);
                    y = BOHR2NM * boost::lexical_cast<double>(split[2]);
                    z = BOHR2NM * boost::lexical_cast<double>(split[3]);
                }
                else if (units == "angstrom") {
                    x = ANGSTROM2NM * boost::lexical_cast<double>(split[1]);
                    y = ANGSTROM2NM * boost::lexical_cast<double>(split[2]);
                    z = ANGSTROM2NM * boost::lexical_cast<double>(split[3]);
                }
                else {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }

                vec pos = vec(x,y,z);

                int rank = boost::lexical_cast<int>(split[5]);

                PolarSite *newPole = new PolarSite(id, name);
                newPole->setRank(rank);
                newPole->setPos(pos);
                poles.push_back(newPole);
                thisPole = newPole;

            }

            // 'P', dipole polarizability
            else if ( split[0] == "P" && split.size() == 2 ) {
                P1 = 1e-3 * boost::lexical_cast<double>(split[1]);
                thisPole->setPs(P1, state);
                useDefaultPs = false;
            }

            // multipole line
            else {

                int lineRank = int( sqrt(thisPole->getQs(state).size()) + 0.5 );

                if (lineRank == 0) {
                    Q0_total += boost::lexical_cast<double>(split[0]);
                }

                for (int i = 0; i < split.size(); i++) {

                    double qXYZ = boost::lexical_cast<double>(split[i]);

                    // Convert e*(a_0)^k to e*(nm)^k where k = rank
                    double BOHR2NM = 0.0529189379;
                    qXYZ *= pow(BOHR2NM, lineRank); // OVERRIDE

                    Qs.push_back(qXYZ);

                }
                thisPole->setQs(Qs, state);
            }

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << filename << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl << "... ... Reading " << filename <<
                    ": Q0(Total) = " << Q0_total << flush;

    if (useDefaultPs) {

        cout << endl << "... ... ... NOTE Using default Thole polarizabilities "
             << "for charge state " << state << ". ";

        vector< PolarSite* > ::iterator pol;
        for (pol = poles.begin(); pol < poles.end(); ++pol) {
            string elem = (*pol)->getName();
            double alpha = 0.0;
            // Original set of Thole polarizabilites
            if      (elem == "C") { alpha = 1.75e-3;  } // <- conversion from
            else if (elem == "H") { alpha = 0.696e-3; } //    A³ to nm³ = 10⁻³
            else if (elem == "N") { alpha = 1.073e-3; }
            else if (elem == "O") { alpha = 0.837e-3; }
            else if (elem == "S") { alpha = 2.926e-3; }
            // Different set of Thole polarizabilities
            //if      (elem == "C") { alpha = 1.334e-3; } // <- conversion from
            //else if (elem == "H") { alpha = 0.496e-3; } //    A³ to nm³ = 10⁻³
            //else if (elem == "N") { alpha = 1.073e-3; }
            //else if (elem == "O") { alpha = 0.837e-3; }
            //else if (elem == "S") { alpha = 3.300e-3; }
            else { throw runtime_error("No polarizability given "
                                       "for polar site type " + elem + ". "); }
            (*pol)->setPs(alpha, state);
        }
    }

    return poles;
}


bool XMP::EvaluateFrame(Topology *top) {

    Collect_XML(_xml_file, top); // <- Allocate polar sites to fragments    
    Collect_EMP(_emp_file, top); // <- Collect segs .mps => Background .mps
    Collect_JOB(_job_file, top); // <- Collect jobs .mps => Foreground .mps
    Collect_MPS(top);            // <- Parse all .mps files (fore-, background)


    // ++++++++++++++++++++ //
    // Ridigidfy + Estatify //
    // ++++++++++++++++++++ //

    // Rigidify segments, fragments
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else {
        cout << endl
             << "... ... System is already rigidified.";
    }

    // Create polar sites
    if (top->isEStatified() == false) {
        this->Create_MPOLS(top);
        cout << endl
             << "... ... Created " << top->PolarSites().size()
             << " multipole sites. "
             << flush;
    }
    else {
        cout << endl
             << "... ... System is already estatified. "
             << "This prohibits running XMultipole. "
             << flush;
    }

    // To check rotation into global frame
    if (_pdb_check != "") {
        cout << endl
             << "... ... Writing polar-site coordinates to "
             << _pdb_check << ". "
             << flush;

        string mpNAME   = _pdb_check;
        FILE *mpPDB     = NULL;
        mpPDB           = fopen(mpNAME.c_str(), "w");

        vector<Segment*>::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {
            (*sit)->WritePDB(mpPDB, "Multipoles", "MD");
        }
        fclose(mpPDB);
    }


}


}}

#endif

