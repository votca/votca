#ifndef XMULTIPOLE2_H
#define XMULTIPOLE2_H


#include <votca/xtp/qmcalculator.h>

namespace votca { namespace xtp {

class XMP : public QMCalculator
{

public:

    XMP() {};
   ~XMP() {};

    string              Identify() { return "XMultipole"; }
    void                Initialize(Property *);
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
    // XJob (distr. among threads) //
    // +++++++++++++++++++++++++++ //

    class XJob
    {
    public:

        XJob(int id,      string tag,  string job_type, int site_id,
             int pairId,  int seg1Id,  int seg2Id,
             string mps1, string mps2, Topology *top) :

             _id(id),         _tag(tag),       _pairId(pairId),
             _type(job_type), _seg1Id(seg1Id), _seg2Id(seg2Id),
             _mps1(mps1),     _mps2(mps2),     _start_from_cpt(false),
             _site_id(site_id) {

             _seg1 = top->getSegment(seg1Id);
             _seg2 = top->getSegment(seg2Id);      

             // JOB TYPE: PAIR
             if (job_type == "pair") {

                 // Watch out for periodic-boundary correction
                 _pair = top->NBList().FindPair(_seg1,_seg2);

                 if (_pair == NULL) {
                    cout << endl
                         << "WARNING: No such pair " << pairId
                         << "(SEG1 " << seg1Id << ", SEG2 " << seg2Id << "). "
                         << flush;
                    throw runtime_error("Pair specs do not match topology.");
                 }

                 assert ( pairId == _pair->getId() );                           // Inconsistency in job input w.r.t. top

                 _center = ( _pair->Seg1PbCopy()->Atoms().size()
                           * _pair->Seg1PbCopy()->getPos()

                           + _pair->Seg2PbCopy()->Atoms().size()
                           * _pair->Seg2PbCopy()->getPos() )

                           / (_pair->Seg1PbCopy()->Atoms().size()
                           +  _pair->Seg2PbCopy()->Atoms().size() );
             }

             // JOB TYPE: SITE
             else if (job_type == "site") {

                 if      (_site_id == seg1Id) { _center = _seg1->getPos();
                                                _site   = _seg1;           }
                 else if (_site_id == seg2Id) { _center = _seg2->getPos();
                                                _site   = _seg2;           }
                 else    { throw runtime_error("ERROR: "
                       "No such site " + boost::lexical_cast<string>(_site_id) +
                       " in job " + boost::lexical_cast<string>(_id) + ". "); }
             }

             // JOB TYPE: UNKNOWN
             else {
                cout << endl 
                     << "ERROR: Job " << _id << ": Invalid job type "
                     << job_type << ": Should be either 'pair' or 'site'"
                     << endl;
             }
       }

       ~XJob() {};

       int      getId()                     { return _id; }
       string   getTag()                    { return _tag; }
       int      getPairId()                 { return _pairId; }
       int      getSeg1Id()                 { return _seg1Id; }
       int      getSeg2Id()                 { return _seg2Id; }
       string   getMPS1()                   { return _mps1; }
       string   getMPS2()                   { return _mps2; }
       int      getSiteId()                 { return _site_id; }
       string   getType()                   { return _type; }

       void     setType(string typ, int id) { _type = typ; _site_id = id; }
       void     setIter(int iter)           { _iter = iter; }
       void     setSizePol(int size)        { _sizePol = size; }
       void     setSizeShell(int size)      { _sizeShell = size; }
       void     setEnergy(double E_T, double E_i_i, 
                          double E_i_I, double E_I_I, double E_i_O,
                          double E_P, double E_U) {
           _E_Tot = E_T;
           _E_Pair_Pair = E_i_i;
           _E_Pair_Sph1 = E_i_I;
           _E_Sph1_Sph1 = E_I_I;
           _E_Pair_Sph2 = E_i_O;

           _E_PERM      = E_P;
           _E_INDU      = E_U;
       }

       void     WriteInfoLine(FILE *out) {

         if (_type == "pair") {
           fprintf(out, "%5d %-20s  E_TOT %+4.7f E_PAIR_PAIR %+4.7f"
                        " E_PAIR_CUT1 %+4.7f E_CUT1_CUT1 %+4.7f E_PAIR_CUT2 "
                        "%+4.7f E_PERM %+4.7f E_INDU %+4.7f ITER %3d"
                        " SPHERE %4d SHELL %4d CENTER %4.7f %4.7f %4.7f \n",
                        _id, _tag.c_str(), _E_Tot, _E_Pair_Pair, _E_Pair_Sph1,
                        _E_Sph1_Sph1, _E_Pair_Sph2, _E_PERM, _E_INDU, _iter,
                        _sizePol, _sizeShell, _center.getX(), _center.getY(),
                        _center.getZ() );
         }
         else if (_type == "site") {
           fprintf(out, "%5d %-20s  E_TOT %+4.7f E_SITE_SITE %+4.7f"
                        " E_SITE_CUT1 %+4.7f E_CUT1_CUT1 %+4.7f E_SITE_CUT2 "
                        "%+4.7f E_PERM %+4.7f E_INDU %+4.7f ITER %3d"
                        " SPHERE %4d SHELL %4d CENTER %4.7f %4.7f %4.7f \n",
                        _id, _tag.c_str(), _E_Tot, _E_Pair_Pair, _E_Pair_Sph1,
                        _E_Sph1_Sph1, _E_Pair_Sph2, _E_PERM, _E_INDU, _iter,
                        _sizePol, _sizeShell, _center.getX(), _center.getY(),
                        _center.getZ() );
         }

       }

       bool     StartFromCPT()  { return _start_from_cpt; }

       Segment *Seg1()          { return _seg1; }
       Segment *Seg2()          { return _seg2; }
       QMPair  *Pair()          { return _pair; }
       vec     &Center()        { return _center; }

    private:

       int          _id;
       string       _tag;
       int          _pairId;
       string       _type;
       int          _seg1Id;
       int          _seg2Id;
       string       _mps1;
       string       _mps2;

       bool         _start_from_cpt;
       int          _site_id; // only relevant if type == "site"

       QMPair      *_pair;
       Segment     *_site;
       Segment     *_seg1;
       Segment     *_seg2;
       vec          _center;


       int          _iter;
       int          _sizePol;
       int          _sizeShell;

       double       _E_Tot;
       double       _E_Pair_Pair;
       double       _E_Pair_Sph1;
       double       _E_Sph1_Sph1;
       double       _E_Pair_Sph2;

       double       _E_PERM;
       double       _E_INDU;

    };

    // +++++++++++++++++++++++++++ //
    // Induction Op    (Subthread) //
    // +++++++++++++++++++++++++++ //

    class JobXMP;

    class InduWorker : public Thread
    {
      public:

        InduWorker(int id, Topology *top, XMP *master, JobXMP *forker)
                    : _id(id), /*_top(top),*/ _master(master), _forker(forker)
                  { _actor = XInteractor(top,_master); };

       ~InduWorker() {};

        void Run(void) {

            if (_switch_energy_induce == 1) {
                while (_forker->NextChunkTodo(_id,1)) {
                    this->InterInduce();
                    _forker->OpenChunks(_chunk1, _chunk2);
                }
            }
            else if (_switch_energy_induce == 0) {
                while (_forker->NextChunkTodo(_id,0)) {
                    this->InterEnergy();
                    _forker->OpenChunks(_chunk1, _chunk2);
                }
            }
            else {
                assert(false);
            }


        }

        void InitSpheres(vector< Segment* > *vsegs_cut1,
                         vector< Segment* > *vsegs_cut2,
                         vector< vector<PolarSite*> > *vvpoles_cut1,
                         vector< vector<PolarSite*> > *vvpoles_cut2) {
            _vsegs_cut1 = vsegs_cut1;
            _vsegs_cut2 = vsegs_cut2;
            _vvpoles_cut1 = vvpoles_cut1;
            _vvpoles_cut2 = vvpoles_cut2;

            _vvpoles = vvpoles_cut1;
        }

        void SetSwitch(int energy_induce) {
            _switch_energy_induce = energy_induce;
            if (energy_induce == 0) {
                _actor.ResetEnergy();
                _E_Pair_Pair = 0.0;
                _E_Pair_Sph1 = 0.0;
                _E_Sph1_Sph1 = 0.0;
            }
            
        }

        void Setc1c2(int c1, int c2) {
            _chunk1 = c1;
            _chunk2 = c2;
        }

        void Setnx12ny12(int nx1, int nx2, int ny1, int ny2) {
            _nx1 = nx1;
            _nx2 = nx2;
            _ny1 = ny1;
            _ny2 = ny2;
        }

        void PrintInfo() {
            printf("%1d : nx1 %2d nx2 %2d ny1 %2d ny2 %2d\n",
                    _id, _nx1, _nx2, _ny1, _ny2);
        }

        void InterInduce() {

            for (int i = _nx1;                     i < _nx2; ++i) {
            for (int j = (i >= _ny1) ? i+1 : _ny1; j < _ny2; ++j) {

                for (pit1 = (*_vvpoles)[i].begin(); pit1 < (*_vvpoles)[i].end(); ++pit1) {
                for (pit2 = (*_vvpoles)[j].begin(); pit2 < (*_vvpoles)[j].end(); ++pit2) {

                    _actor.FieldIndu(*(*pit1), *(*pit2));
                }}
            }}
        }

        double       GetEPairPair() { return _E_Pair_Pair; }
        double       GetEPairSph1() { return _E_Pair_Sph1; }
        double       GetESph1Sph1() { return _E_Sph1_Sph1; }
        XInteractor &GetActor()     { return _actor; }

        void InterEnergy() {

            for (int i = _nx1;                     i < _nx2; ++i) {
            for (int j = (i >= _ny1) ? i+1 : _ny1; j < _ny2; ++j) {

                // Site-non-site interaction
                if (this->_forker->_job->getSiteId() == (*_vsegs_cut1)[i]->getId()
                 || this->_forker->_job->getSiteId() == (*_vsegs_cut1)[j]->getId()) {

                    for (pit1 = (*_vvpoles_cut1)[i].begin();
                         pit1 < (*_vvpoles_cut1)[i].end();
                         ++pit1) {
                    for (pit2 = (*_vvpoles_cut1)[j].begin();
                         pit2 < (*_vvpoles_cut1)[j].end();
                         ++pit2) {

                        _E_Pair_Sph1 += _actor.EnergyInter(*(*pit1),*(*pit2));
                    }}
                }

                // Non-site-non-site interaction
                else {
                    for (pit1 = (*_vvpoles_cut1)[i].begin();
                         pit1 < (*_vvpoles_cut1)[i].end();
                         ++pit1) {
                    for (pit2 = (*_vvpoles_cut1)[j].begin();
                         pit2 < (*_vvpoles_cut1)[j].end();
                         ++pit2) {

                        _E_Sph1_Sph1 += _actor.EnergyInter(*(*pit1),*(*pit2));
                    }}
                }
            }}
        }


      private:

          int                                   _id;
          //Topology                             *_top;
          XMP                                  *_master;
          JobXMP                               *_forker;
          XInteractor                           _actor;
          vector< vector<PolarSite*> >         *_vvpoles;

          vector< Segment* >                   *_vsegs_cut1;
          vector< Segment* >                   *_vsegs_cut2;
          vector< vector<PolarSite*> >         *_vvpoles_cut1;
          vector< vector<PolarSite*> >         *_vvpoles_cut2;

          int _nx1, _nx2;
          int _ny1, _ny2;
          int _chunk1, _chunk2;

          int _switch_energy_induce;

          vector< Segment* >              ::iterator      seg1;
          vector< Segment* >              ::iterator      seg2;
          vector< vector<PolarSite*> >    ::iterator      sit1;
          vector< vector<PolarSite*> >    ::iterator      sit2;
          vector< PolarSite* >            ::iterator      pit1;
          vector< PolarSite* >            ::iterator      pit2;

          double _E_Pair_Pair;
          double _E_Pair_Sph1;
          double _E_Sph1_Sph1;

    };

    // +++++++++++++++++++++++++++ //
    // Job Operator (Thread class) //
    // +++++++++++++++++++++++++++ //

    class JobXMP : public Thread
    {
    public:

        JobXMP(int id,     Topology *top,        XMP *master)
                : _id(id),          /*_top(top),*/       _master(master)
                   { _actor = XInteractor(top,_master); };

       ~JobXMP() {};

        int         getId() { return _id; }
        void        setId(int id) { _id = id; }

        void        InitSlotData(Topology *top);
        void        Run(void);
        void        EvalJob(Topology *top, XJob *job);

        int         Induce(int state, XJob *job);
        double      Energy(int state, XJob *job);
        double      EnergyStatic(int state, XJob *job);


        void        ClearTodoTable() {
            for (unsigned int i = 0; i < _xy_done.size(); ++i) {
            for (unsigned int j = 0; j < _xy_done[i].size(); ++j) {
                    _xy_done[i][j] = false;
            }}
        }

        void        InitChunks() {

            _nx1.clear();
            _nx2.clear();
            _ny1.clear();
            _ny2.clear();

            _xy_done.clear();
            _chunks_avail.clear();


            int T = this->_master->_subthreads;         // Threads
            int C = T * 2;                              // Chunks
            int N = _polsPolSphere.size();              // Elements
            int nr = N % C;                             // Rest size
            int nt = (N-nr) / C;                        // Chunk size

            assert (N == C*nt + nr);

            for (int id = 0; id < C+1; ++id) {
                _nx1.push_back( vector<int>(C+1,0) );
                _nx2.push_back( vector<int>(C+1,0) );
                _ny1.push_back( vector<int>(C+1,0) );
                _ny2.push_back( vector<int>(C+1,0) );

                _xy_done.push_back( vector<bool>(C+1,false) );
                _chunks_avail.push_back(true);
            }

            for (int col = 0; col < C+1; ++col) {
                for (int row = 0; row < C+1; ++row) {

                    if (col < row) {
                        _nx1[row][col] = 0; _ny1[row][col] = 0;
                        _nx2[row][col] = 0; _ny2[row][col] = 0;
                    }
                    else if (col == C && row == C) {
                        _nx1[row][col] = C*nt;
                        _nx2[row][col] = C*nt + nr;
                        _ny1[row][col] = C*nt;
                        _ny2[row][col] = C*nt + nr;
                    }
                    else if (col == C && row < C) {
                        _nx1[row][col] = row*nt;
                        _nx2[row][col] = (row+1)*nt;
                        _ny1[row][col] = C*nt;
                        _ny2[row][col] = C*nt + nr;
                    }
                    else {
                        _nx1[row][col] = row*nt;
                        _nx2[row][col] = (row+1)*nt;
                        _ny1[row][col] = col*nt;
                        _ny2[row][col] = (col+1)*nt;
                    }
                }
            }

            if (T > 1) {
                this->_master->_coutMutex.Lock();
                printf("\n\nTHREAD %1d MESH LOAD: "
                       "NST%1d C%1d N%1d nt%1d nr%1d\n", _id, T, C, N, nt, nr);
                for (int id = 0; id < C+1; ++id) {
                    for (int run = 0; run < C+1; ++run) {
                        printf("--------+");
                    }
                    printf("\n");
                    for (int run = 0; run < C+1; ++run) {
                        printf("%3d %3d |", _nx1[id][run], _ny1[id][run]);
                    }
                    printf("\n");
                    for (int run = 0; run < C+1; ++run) {
                        printf("%3d %3d |", _nx2[id][run], _ny2[id][run]);
                    }
                    printf("\n");
                }
                for (int run = 0; run < C+1; ++run) {
                    printf("--------+");
                }
                printf("\n");
                this->_master->_coutMutex.Unlock();
            }
        }


        void        OpenChunks(int c1, int c2) {
            _chunks_avail[c1] = true;
            _chunks_avail[c2] = true;            
        }

        bool        NextChunkTodo(int indu_id, int switch_energy_induce) {

            _alloc_chunk.Lock();

            bool todo = false;

            while (true) {

                for (unsigned int i = 0; i < _xy_done.size(); ++i) {
                for (unsigned int j = i; j < _xy_done[i].size(); ++j) {
                    if (!_xy_done[i][j]) {
                        todo = true;
                        if (_chunks_avail[i] && _chunks_avail[j]) {

                            if (switch_energy_induce == 1) {
                                _chunks_avail[i] = false;
                                _chunks_avail[j] = false;
                            }

                            _xy_done[i][j] = true;
                            _indus[indu_id]->Setc1c2(i,j);
                            _indus[indu_id]->Setnx12ny12(_nx1[i][j],
                                                         _nx2[i][j],
                                                         _ny1[i][j],
                                                         _ny2[i][j]);

                            _alloc_chunk.Unlock();
                            return todo;
                        }
                        else { ; }
                    }
                }}

                if (!todo) { break; }

            }

            _alloc_chunk.Unlock();
            return todo;
        }
        



    public:

        int                           _id;
        Topology                     *_top;
        XMP                          *_master;
        XJob                         *_job;

        vector< Segment* >           _segsPolSphere;  // Segments    in c/o 0-1
        vector< Segment* >           _segsOutSphere;  // Segments    in c/0 1-2
        vector< vector<PolarSite*> > _polsPolSphere;  // Polar sites in c/o 0-1
        vector< vector<PolarSite*> > _polsOutSphere;  // Polar sites in c/o 1-2
        vector< vector<PolarSite*> > _polarSites;     // Copy of top polar sites
        vector< vector<PolarSite*> > _polarSites_job; // Adapted to job specs
        XInteractor                  _actor;

        // Manage induction workers
        vector< InduWorker* >        _indus;

        Mutex                        _alloc_chunk;
        vector< bool >               _chunks_avail;

        vector< vector<bool> >       _xy_done;
        
        vector< vector<int> >        _nx1;
        vector< vector<int> >        _nx2;
        vector< vector<int> >        _ny1;
        vector< vector<int> >        _ny2;


    };

    XJob               *RequestNextJob(int id, Topology *top);
    void                LockCout() { _coutMutex.Lock(); }
    void                UnlockCout() { _coutMutex.Unlock(); }





private:

    // ++++++++++++++++++++++++++++++++++++++++ //
    // MULTIPOLE ALLOCATION TO SEGMENTS / PAIRS //
    // ++++++++++++++++++++++++++++++++++++++++ //

    string                         _job_file;
    string                         _emp_file;
    string                         _xml_file;

    // Job info : JOB_ID PAIR_ID MPS_1 MPS_2 TAG
    vector<XJob*>                  _XJobs;
    vector<XJob*>::iterator        _nextXJob;
    Mutex                          _nextJobMutex;
    Mutex                          _coutMutex;
    Mutex                          _logMutex;
    bool                           _maverick;

    vector<int>                    _jobIds;
    map<int,string>                _jobId_jobTag;
    map<int,int>                   _jobId_pairId;
    map<int,string>                _jobId_mpsFile1;
    map<int,string>                _jobId_mpsFile2;
    map<int, pair<int,int> >       _pairId_seg1Id_seg2Id;
    vector<int>::iterator          _nextJob;

    // Allocate .mps files to segments (n 0, e -1, h +1)
    map<int,string>                 _segId_mpsFile_n;
    map<int,string>                 _segId_mpsFile_e;
    map<int,string>                 _segId_mpsFile_h;

    // Store polar site containers, one for each .mps file
    map<string,vector<PolarSite*> > _mpsFile_pSites;
    map<string,vector<PolarSite*> > _mpsFile_pSites_job;

    // Allocate polar sites to fragments in segments
    map<string, bool>               _map2md;
    map<string, vector<int> >       _alloc_frag_mpoleIdx;
    map<string, vector<string> >    _alloc_frag_mpoleName;
    map<string, vector<int> >       _alloc_frag_trihedron;
    map<string, vector<double> >    _alloc_frag_weights;
    //map<string, vector<int> >       _alloc_frag_trihedron_mps;
    //map<string, vector<double> >    _alloc_frag_weights_mps;

    string                          _pdb_check;
    bool                            _write_chk;
    string                          _write_chk_suffix;
    bool                            _chk_split_dpl;
    double                          _chk_dpl_spacing;
    string                          _chk_format;


    // ++++++++++++++++++++++++++++++++++++++++ //
    // INDUCTION + ENERGY EVALUATION            //
    // ++++++++++++++++++++++++++++++++++++++++ //

    // Control options
    bool                            _induce;
    bool                            _induce_intra_pair;
    int                             _subthreads;

    // Interaction parameters
    bool                            _useCutoff;
    double                          _cutoff;
    double                          _cutoff2;
    bool                            _useExp;
    double                          _aDamp;
    bool                            _useScaling;
    vector<double>                  _scale1;

    // Convergence parameters
    float                           _wSOR_N;
    float                           _wSOR_C;
    double                          _epsTol;
    int                             _maxIter;

    // Logging
    string                          _outFile;
    bool                            _energies2File;
    map<int,vector<int> >           _log_seg_iter;
    map<int,int>                    _log_seg_sphereSize;
    map<int,vec>                    _log_seg_com;


};

// ========================================================================== //
//                        XMULTIPOLE MEMBER FUNCTIONS                         //
// ========================================================================== //


void XMP::Initialize(Property *opt) {

    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;
    

    string key = "options.xmultipole.multipoles";

        if ( opt->exists(key) ) {
            _xml_file = opt->get(key).as< string >();
        }

    key = "options.xmultipole.control";

        if ( opt->exists(key+".job_file")) {
            _job_file   = opt->get(key+".job_file").as<string>();
        }
        else {
            _job_file   = opt->get(key+".job_file").as<string>();
        }

        if ( opt->exists(key+".emp_file")) {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }
        else {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }

        if ( opt->exists(key+".output") ) {
            _outFile = opt->get(key+".output").as< string >();
            _energies2File = true;
        }
        else { _energies2File = false; }

        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<string>();
        }
        else { _pdb_check = ""; }

        if (opt->exists(key+".write_chk")) {
            _write_chk_suffix = opt->get(key+".write_chk").as<string>();
            _write_chk = true;
        }
        else { _write_chk = false; }

        if (opt->exists(key+".format_chk")) {
            _chk_format = opt->get(key+".format_chk").as<string>();
        }
        else { _chk_format = "xyz"; }

        if (opt->exists(key+".split_dpl")) {
            _chk_split_dpl = (opt->get(key+".split_dpl").as<int>() == 1) ?
                         true : false;
        }
        else { _chk_split_dpl = true; }

        if (opt->exists(key+".dpl_spacing")) {
            _chk_dpl_spacing = opt->get(key+".dpl_spacing").as<double>();
        }
        else {
            _chk_dpl_spacing = 1.0e-6;
        }


    key = "options.xmultipole.tholemodel";

        if ( opt->exists(key+".induce") ) {
            int induce = opt->get(key+".induce").as< int >();
            _induce = (induce == 0) ? false : true;
        }
        else { _induce = true; }

        if ( opt->exists(key+".induce_intra_pair") ) {
            int induce = opt->get(key+".induce_intra_pair").as< int >();
            _induce_intra_pair = (induce == 0) ? false : true;
        }
        else { _induce_intra_pair = true; }

        if ( opt->exists(key+".cutoff1") ) {
            _cutoff = opt->get(key+".cutoff1").as< double >();
            if (_cutoff) { _useCutoff = true; }
        }
        if ( opt->exists(key+".cutoff2") ) {
            _cutoff2 = opt->get(key+".cutoff2").as< double >();
        }
        else {
            _cutoff2 = _cutoff;
        }
        if ( opt->exists(key+".subthreads") ) {
            _subthreads = opt->get(key+".subthreads").as< double >();
        }
        else {
            _subthreads = 1;
        }
        if ( opt->exists(key+".exp_damp") ) {
            _aDamp = opt->get(key+".exp_damp").as< double >();
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

    key = "options.xmultipole.convergence";

        if ( opt->exists(key+".wSOR_N") ) {
            _wSOR_N = opt->get(key+".wSOR_N").as< float >();
        }
        else { _wSOR_N = 0.75; }
        if ( opt->exists(key+".wSOR_C") ) {
            _wSOR_C = opt->get(key+".wSOR_C").as< float >();
        }
        else { _wSOR_C = 0.75; }

        if ( opt->exists(key+".max_iter") ) {
            _maxIter = opt->get(key+".max_iter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; }    

}


void XMP::Collect_XML(string xml_file, Topology *top) {

    cout << endl 
         << "... ... ... Allocate polar sites to fragments. "
         << flush;

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

                // Local frame for polar sites
                vector<int> trihedron_mps;
                if ((*fragit)->exists("localframe_mps")) {
                   cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Defining distinct local frame for polar sites."
                         << flush;
                   trihedron_mps = (*fragit)->get("localframe_mps")
                                         .as< vector<int> >();
                }
                else {
                   trihedron_mps = (*fragit)->get("localframe")
                                         .as< vector<int> >();
                }
                
                // Mapping weights for polar sites
                vector<double> weights_mps;
                if ((*fragit)->exists("weights_mps")) {
                    cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Using distinct weights for polar sites."
                         << flush;
                   weights_mps = (*fragit)->get("weights_mps")
                                       .as< vector<double> >();
                }
                else {
                   weights_mps = (*fragit)->get("weights")
                                       .as< vector<double> >();
                }

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
                //alloc_frag_trihedron_mps[mapKeyName] = trihedron_mps;
                //alloc_frag_weights_mps[mapKeyName] = weights_mps;
                alloc_frag_trihedron[mapKeyName] = trihedron_mps;
                alloc_frag_weights[mapKeyName] = weights_mps;
            }
        }
    }

    _alloc_frag_mpoleIdx        = alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName       = alloc_frag_mpoleName;
    //_alloc_frag_trihedron_mps   = alloc_frag_trihedron_mps;
    //_alloc_frag_weights_mps     = alloc_frag_weights_mps;
    _alloc_frag_trihedron       = alloc_frag_trihedron;
    _alloc_frag_weights         = alloc_frag_weights;
}


void XMP::Collect_JOB(string job_file, Topology *top) {

    //QMNBList &nblist = top->NBList();

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
// # JOB_ID TAG  PAIR_ID SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS  (TYPE  SITE)
//   1      E_CT 3819    182     C60       c60.mps  392     DCV       dcv.mps   (site  392)

            int jobId       = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            int pairId      = boost::lexical_cast<int>(split[2]);

            int seg1Id      = boost::lexical_cast<int>(split[3]);
            string seg1Name = split[4];
            string seg1mps  = split[5];

            int seg2Id      = boost::lexical_cast<int>(split[6]);
            string seg2Name = split[7];
            string seg2mps  = split[8];

            string job_type = "pair";
            int    energy_site_id = -1;
            if (split.size() == 11) {
                job_type = split[9];
                energy_site_id = boost::lexical_cast<int>(split[10]);
            }

            //Segment *seg1   = top->getSegment(seg1Id);
            //Segment *seg2   = top->getSegment(seg2Id);

            _XJobs.push_back(new XJob(jobId,  tag,    job_type, energy_site_id,
                                      pairId, seg1Id, seg2Id,   seg1mps,
                                      seg2mps, top));

            _XJobs.back()->setType(job_type, energy_site_id);

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl 
         << "... ... ... Registered " << _XJobs.size() << " jobs. "
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


    // Job Seg1 Seg2
    vector<XJob*> :: iterator jit;
    for (jit = _XJobs.begin();
         jit < _XJobs.end();
         ++jit) {

        string mpsfile1 = (*jit)->getMPS1();
        string mpsfile2 = (*jit)->getMPS2();

        if (_mpsFile_pSites_job.count(mpsfile1) > 0 ) { ; }
        else { _mpsFile_pSites_job[mpsfile1] = Parse_GDMA(mpsfile1, 0); }

        if (_mpsFile_pSites_job.count(mpsfile2) > 0 ) { ; }
        else { _mpsFile_pSites_job[mpsfile2] = Parse_GDMA(mpsfile2, 0); }

    }

//    map<string, vector<PolarSite*> > ::iterator it;
//    for (it = _mpsFile_pSites_job.begin();
//         it != _mpsFile_pSites_job.end();
//         ++it) {
//         cout << endl << "KEY_" << it->first << flush;
//    }


    cout << endl
         << "... ... ... Parsed " << _mpsFile_pSites.size()
         << " .mps files from " << _emp_file
         << ", " << _mpsFile_pSites_job.size()
         << " .mps files from " << _job_file
         << flush;
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
    
    // Warning: Direct mapping of k > 0 multipoles to MD coordinates
    bool print_huge_map2md_warning = false;
    
    // Log warning: Symmetry = 1 and k > 0 multipoles.
    map<string,bool> warned_symm_idkey;  
    
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

        for (unsigned int i = 0; i < pols_n.size(); i++) {
            pols_n[i]->setQs( pols_h[i]->getQs(+1), +1 );
            pols_n[i]->setPs( pols_h[i]->getPs(+1), +1 );
        }
        for (unsigned int i = 0; i < pols_n.size(); ++i) {
            pols_n[i]->setQs(pols_e[i]->getQs(-1), -1 );
            pols_n[i]->setPs(pols_e[i]->getPs(-1), -1 );
        }

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            // Retrieve polar-site data for this fragment
            string idkey                     = frag->getName() + seg->getName()
                                             + seg->getMolecule()->getName();
            vector<int> polesInFrag          = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag       = _alloc_frag_mpoleName.at(idkey);
            vector<double> weightsInFrag     = _alloc_frag_weights.at(idkey);

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
                }
                
                trihedron_ints  = frag->getTrihedron();
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    vector< Atom* > ::iterator ait;
                    for (ait = frag->Atoms().begin();
                         ait < frag->Atoms().end();
                         ++ait) {
                        if ((*ait)->getQMId() == (*iit)) {
                            trihedron_atm.push_back(*ait);
                        }
                    }
                }


                int symmetry = trihedron_pol.size();
                assert (trihedron_pol.size() <= trihedron_atm.size() );                

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

                else if (symmetry == 2) {

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

                else if (symmetry == 1) {
                    
                    if (!warned_symm_idkey[idkey]) {
                        cout << endl << "... ... ... "
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << flush;
                        warned_symm_idkey[idkey] = true;
                    }

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
                for (unsigned int i = 0; i < polesInFrag.size(); ++i) {

                    double weight = weightsInFrag[i];
                    
                    vec pos = pols_n[polesInFrag[i]-1]->getPos();
                    
                    CoMP += weight*pos;
                    W += weight;

                }
                CoMP /= W;

                translateMP2MD = frag->getCoMD() - CoMP;

            }            

            // Create polar sites 
            for (unsigned int i = 0; i < polesInFrag.size(); i++) {

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
                        print_huge_map2md_warning = true;
                    }
                }
            }
        }
    }

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

    top->setIsEStatified(true);

}


vector<PolarSite*> XMP::Map_MPols_To_Seg(vector<PolarSite*> &pols_n, Segment *seg) {

    bool print_huge_map2md_warning = false;

    vector<PolarSite*> return_pols;
    return_pols.reserve(pols_n.size());

    //int segId                   = seg->getId();
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
            }

            trihedron_ints  = frag->getTrihedron();
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                vector< Atom* > ::iterator ait;
                for (ait = frag->Atoms().begin();
                     ait < frag->Atoms().end();
                     ++ait) {
                    if ((*ait)->getQMId() == (*iit)) {
                        trihedron_atm.push_back(*ait);
                    }
                }
            }

            int symmetry = trihedron_pol.size();
            assert (trihedron_pol.size() <= trihedron_atm.size() );

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

            else if (symmetry == 2) {

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

            else if (symmetry == 1) {

                //cout << endl
                //     << "WARNING: Symmetry = 1 for fragment "
                //     << frag->getName() << ": This will generate artifacts "
                //     << "when mapping higher-rank multipoles (dipoles, ..)."
                //     << endl;

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
            for (unsigned int i = 0; i < polesInFrag.size(); ++i) {

                double weight = weightsInFrag[i];

                vec pos = pols_n[polesInFrag[i]-1]->getPos();

                CoMP += weight*pos;
                W += weight;

            }
            CoMP /= W;

            translateMP2MD = frag->getCoMD() - CoMP;

        }

        // Create polar sites
        for (unsigned int i = 0; i < polesInFrag.size(); i++) {

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
                    print_huge_map2md_warning = true;
                }
            }

            newSite->Charge(0);
            
            return_pols.push_back(newSite);

        }
    } // End loop over fragments

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

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

                for (unsigned int i = 0; i < split.size(); i++) {

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


    printf("\n... ... ... Reading %-25s -> Q0(Sum) = %+1.3f ",
                          filename.c_str(),       Q0_total);

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

    cout << endl
         << "... ... Load multipole definition, collect jobs: "
         << flush;

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

    
    // +++++++++++++++++++++++++++++++++ //
    // Create + start threads (Job Ops)  //
    // +++++++++++++++++++++++++++++++++ //

    // Convert threads into subthreads if beneficial
    if (_XJobs.size() < _nThreads) {
        _subthreads = (_nThreads - _XJobs.size()) / _XJobs.size() + 1;
        _nThreads   = _XJobs.size();

        cout << endl << "... ... "
             << "Converted threads into subthreads to increase efficiency: "
             << "NT = " << _nThreads << ", NST = " << _subthreads
             << flush;
    }

    vector<JobXMP*> jobOps;    

    _nextXJob = _XJobs.begin();

    for (unsigned int id = 0; id < _nThreads; id++) {
        JobXMP *newOp = new JobXMP(id, top, this);
        jobOps.push_back(newOp);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        jobOps[id]->InitSlotData(top);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        jobOps[id]->Start();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        jobOps[id]->WaitDone();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        delete jobOps[id];
    }

    jobOps.clear();


    FILE *out;
    out = fopen(this->_outFile.c_str(), "w");
    vector<XJob*> :: iterator jit;
    for (jit = _XJobs.begin(); jit < _XJobs.end(); ++jit) {
        (*jit)->WriteInfoLine(out);
    }
    fclose(out);

    return true;
    
}



// ========================================================================== //
//                            JOBXMP MEMBER FUNCTIONS                         //
// ========================================================================== //


void XMP::JobXMP::InitSlotData(Topology *top) {

    vector< Segment* > ::iterator sitRef;
    vector< vector<PolarSite*> > ::iterator sitNew;
    vector< PolarSite* > ::iterator pitRef;
    vector< PolarSite* > ::iterator pitNew;

    _polarSites.resize(top->Segments().size());
    assert(top->Segments().size() == _polarSites.size());

    for (sitRef = top->Segments().begin(), sitNew = _polarSites.begin();
         sitRef < top->Segments().end();
         ++sitRef, ++sitNew) {

        (*sitNew).resize((*sitRef)->PolarSites().size());

        for (pitRef = (*sitRef)->PolarSites().begin(),
             pitNew = (*sitNew).begin();
             pitRef < (*sitRef)->PolarSites().end();
             ++pitRef, ++ pitNew) {

            *pitNew = new PolarSite();
            (*pitNew)->ImportFrom(*pitRef, "full");
            (*pitNew)->Charge(0);
        }
    }
}


void XMP::JobXMP::Run(void) {

    while (true) {
        _job = _master->RequestNextJob(_id, _top);

        if (_job == NULL) { break; }
        else { this->EvalJob(_top, _job); }
    }
}


void XMP::JobXMP::EvalJob(Topology *top, XJob *job) {

    //double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    // ++++++++++++++++++++++++++ //
    // Adapt polar sites          //
    // ++++++++++++++++++++++++++ //
    _polarSites_job = _polarSites;

    int subs_here1 = job->getSeg1Id();
    int subs_here2 = job->getSeg2Id();

    vector<PolarSite*> subs1_raw = _master->_mpsFile_pSites_job[job->getMPS1()];
    vector<PolarSite*> subs2_raw = _master->_mpsFile_pSites_job[job->getMPS2()];
    vector<PolarSite*> subs1 = _master->Map_MPols_To_Seg(subs1_raw,job->Seg1());
    vector<PolarSite*> subs2 = _master->Map_MPols_To_Seg(subs2_raw,job->Seg2());

    _polarSites_job[subs_here1-1] = subs1;
    _polarSites_job[subs_here2-1] = subs2;
    
    

    // ++++++++++++++++++++++++++ //
    // Define polarization sphere //
    // ++++++++++++++++++++++++++ //

    vec center = job->Center();
    // vec center = job->Seg1()->getPos(); // UNCOMMENT TO COMPARE TO EMULTIPOLE

    this->_segsPolSphere.clear(); // <- Segments    within cutoff
    this->_segsOutSphere.clear(); // <- Segments    within cutoff1, cutoff2
    this->_polsPolSphere.clear(); // <- Polar sites within cutoff
    this->_polsOutSphere.clear(); // <- Polar sites within cutoff1, cutoff2

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

        double r12 = abs(_top->PbShortestConnect((*sit)->getPos(), center));

        // Always add pair-job segments to polSphere, even for cut-off = 0.0
        if ( (*sit)->getId() == job->getSeg1Id()
          || (*sit)->getId() == job->getSeg2Id()) {
            if   (job->getType() == "pair") { r12 = -1; }
            else                            { ; }
        }

        if      ( r12 > _master->_cutoff2) { continue; }

        else if ( r12 > _master->_cutoff ) {
            _segsOutSphere.push_back(*sit);
            _polsOutSphere.push_back( _polarSites_job[(*sit)->getId() - 1] );
        }
        else {
            _segsPolSphere.push_back(*sit);
            _polsPolSphere.push_back( _polarSites_job[(*sit)->getId() - 1] );
        }
    }

    if (_master->_maverick) {
        cout << endl
             << "... ... ... Segments in polarization sphere: "
             << _segsPolSphere.size()
             << "; segments in static shell: "
             << _segsOutSphere.size()
             << flush;
    }

    job->setSizePol(_polsPolSphere.size());
    job->setSizeShell(_polsOutSphere.size());

//    FILE *out;
//    string shellFile = "OuterShell.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    for (sit = _segsOutSphere.begin(); sit < _segsOutSphere.end(); ++sit) {
//        (*sit)->WritePDB(out, "Multipoles", "");
//    }
//    fclose(out);
//
//    shellFile = "InnerShell.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    for (sit = _segsPolSphere.begin(); sit < _segsPolSphere.end(); ++sit) {
//        (*sit)->WritePDB(out, "Multipoles", "");
//    }
//    fclose(out);
//
//    shellFile = "Pair.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    job->Seg1()->WritePDB(out, "Multipoles", "");
//    job->Seg2()->WritePDB(out, "Multipoles", "");
//    fclose(out);


    // ++++++++++++++++++++++++++ //
    // (De-)polarize, charge to N //
    // ++++++++++++++++++++++++++ //

    if (job->StartFromCPT()) {

        if (_master->_maverick) {
            cout << endl
                 << "... ... ... Loading induced dipoles from .cpt file. "
                 << flush;
        }
        assert(false); // Load induced dipole moments from file
    }

    else {

        vector< vector<PolarSite*> > ::iterator sit;
        vector< PolarSite* >         ::iterator pit;

        // Depolarize inner sphere
        for (sit = _polsPolSphere.begin(); sit < _polsPolSphere.end(); ++sit) {
        for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }}

        // Depolarize outer shell
        for (sit = _polsOutSphere.begin(); sit < _polsOutSphere.end(); ++sit) {
        for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }}
    }

    // +++++++++++++++++ //
    // Induction workers //
    // +++++++++++++++++ //

    for (int id = 0; id < this->_master->_subthreads; ++id) {
        InduWorker *newIndu = new InduWorker(id,this->_top,this->_master,this);
        _indus.push_back(newIndu);
        newIndu->InitSpheres(&_segsPolSphere,&_segsOutSphere,
                             &_polsPolSphere,&_polsOutSphere);
        newIndu->SetSwitch(1);
    }

    this->InitChunks();

    // ++++++++++++++++++++++++++ //
    // Compute state energy       //
    // ++++++++++++++++++++++++++ //

    //double  E_state  = 0.0;
    int     iter     = 0;
    int     state    = 0;

    if (_master->_induce) iter      = this->Induce(state, job);
    if (_master->_induce) /*E_state   =*/ (void)this->Energy(state, job);
    else                  /*E_state   =*/ (void)this->EnergyStatic(state, job);

    job->setIter(iter);



    // ++++++++++++++++++++++++++ //
    // Write Checkpoint File      //
    // ++++++++++++++++++++++++++ //

    if (_master->_write_chk) {

        string dotsuffix = "";

        if (_master->_chk_format == "gaussian") {
            dotsuffix = ".com";
        }
        else if (_master->_chk_format == "xyz") {
            dotsuffix = ".xyz";
        }

        FILE *out;
        string chk_file = job->getTag()+_master->_write_chk_suffix+dotsuffix;
        out = fopen(chk_file.c_str(), "w");

        vector< PolarSite* >         ::iterator pit;
        int pcount = 0;

        // Save coordinates of central pair
        if (job->getType() == "pair") {
            for (sit = _segsPolSphere.begin(); sit < _segsPolSphere.end(); ++sit) {

                int segId = (*sit)->getId();

                if (segId == job->getSeg1Id() || segId == job->getSeg2Id()) {

                    vec pb_shift = job->Center() - (*sit)->getPos()
                         - top->PbShortestConnect((*sit)->getPos(), job->Center());

                    for (pit = _polarSites_job[segId-1].begin();
                         pit < _polarSites_job[segId-1].end();
                         ++pit, ++pcount) {
                        (*pit)->WriteXyzLine(out, pb_shift, _master->_chk_format);
                    }
                }
            }
        }
        else {

            if (job->getSiteId() == job->getSeg2Id()) {
                // Reverse order
                for (sit = _segsPolSphere.end() - 1; sit > _segsPolSphere.begin(); --sit) {

                    int segId = (*sit)->getId();

                    if (segId == job->getSeg1Id() || segId == job->getSeg2Id()) {

                        vec pb_shift = job->Center() - (*sit)->getPos()
                             - top->PbShortestConnect((*sit)->getPos(), job->Center());

                        for (pit = _polarSites_job[segId-1].begin();
                             pit < _polarSites_job[segId-1].end();
                             ++pit, ++pcount) {
                            if (segId == job->getSiteId()) {
                                (*pit)->WriteXyzLine(out, pb_shift, _master->_chk_format);
                            }
                            else {
                                (*pit)->WriteChkLine(out, pb_shift,
                                               _master->_chk_split_dpl,
                                               _master->_chk_format,
                                               _master->_chk_dpl_spacing);
                            }
                        }
                    }
                }
            }
            else if (job->getSiteId() == job->getSeg1Id()) {
                for (sit = _segsPolSphere.begin(); sit < _segsPolSphere.end(); ++sit) {

                    int segId = (*sit)->getId();

                    if (segId == job->getSeg1Id() || segId == job->getSeg2Id()) {

                        vec pb_shift = job->Center() - (*sit)->getPos()
                             - top->PbShortestConnect((*sit)->getPos(), job->Center());

                        for (pit = _polarSites_job[segId-1].begin();
                             pit < _polarSites_job[segId-1].end();
                             ++pit, ++pcount) {
                            if (segId == job->getSiteId()) {
                                (*pit)->WriteXyzLine(out, pb_shift, _master->_chk_format);
                            }
                            else {
                                (*pit)->WriteChkLine(out, pb_shift,
                                               _master->_chk_split_dpl,
                                               _master->_chk_format,
                                               _master->_chk_dpl_spacing);
                            }
                        }
                    }
                }
            }
            else {
                printf("\nERROR: No such segment ID %1d in job ID %1d.",
                        job->getSiteId(), job->getId());
                throw std::runtime_error("Redo job input.");
            }


        }


        if (_master->_chk_format == "gaussian" && job->getType() == "pair") {
            fprintf(out, "\n");
        }

        // Save induction state of polarization sphere
        for (sit = _segsPolSphere.begin(); sit < _segsPolSphere.end(); ++sit) {

            int segId = (*sit)->getId();

            if (segId == job->getSeg1Id() || segId == job->getSeg2Id()) {
                continue;
            }

            vec pb_shift = job->Center() - (*sit)->getPos()
                     - top->PbShortestConnect((*sit)->getPos(), job->Center());


            for (pit = _polarSites_job[segId-1].begin();
                 pit < _polarSites_job[segId-1].end();
                 ++pit, ++pcount) {
                 (*pit)->WriteChkLine(out, pb_shift,
                                           _master->_chk_split_dpl,
                                           _master->_chk_format,
                                           _master->_chk_dpl_spacing);
            }
        }

        // Write point charges of outer sphere
        for (sit = _segsOutSphere.begin(); sit < _segsOutSphere.end(); ++sit) {

            int segId = (*sit)->getId();

            if (segId == job->getSeg1Id() || segId == job->getSeg2Id()) {
                if (job->getType() == "pair") {
                    assert(false);                                              // Central pair in outer shell? No!
                }
                else if (job->getType() == "site") {
                    assert(segId != job->getSiteId());                          // Central site in outer shell? No!
                }                
            }

            vec pb_shift = job->Center() - (*sit)->getPos()
                     - top->PbShortestConnect((*sit)->getPos(), job->Center());

            for (pit = _polarSites_job[segId-1].begin();
                 pit < _polarSites_job[segId-1].end();
                 ++pit, ++pcount) {
                 (*pit)->WriteChkLine(out, pb_shift,
                                           false,
                                           _master->_chk_format,
                                           _master->_chk_dpl_spacing);
            }
        }

        fclose(out);
    }

    // ++++++++++++++++++++++++++ //
    // Clean up polar sites       //
    // ++++++++++++++++++++++++++ //

    vector< PolarSite* > ::iterator cleanit;
    for (cleanit = subs1.begin(); cleanit < subs1.end(); ++cleanit) {
        delete *cleanit;
    }
    for (cleanit = subs2.begin(); cleanit < subs2.end(); ++cleanit) {
        delete *cleanit;
    }
    subs1.clear();
    subs2.clear();

}


int XMP::JobXMP::Induce(int state, XJob *job) {

    for (int id = 0; id < _master->_subthreads; ++id) {
        _indus[id]->SetSwitch(1);
    }

    double wSOR = (state == 0) ? _master->_wSOR_N : _master->_wSOR_C;
    double eTOL = this->_master->_epsTol;
    int    maxI = this->_master->_maxIter;

    // Intra-pair induction ...
    bool   induce_intra_pair = this->_master->_induce_intra_pair;
    // ... change this for jobs of type "site":
    if (job->getType() == "site") { induce_intra_pair = true; }

    vector< vector<PolarSite*> > ::iterator sit1;
    vector< vector<PolarSite*> > ::iterator sit2;
    vector< PolarSite* > ::iterator pit1;
    vector< PolarSite* > ::iterator pit2;
    vector< Segment* > ::iterator seg1;
    vector< Segment* > ::iterator seg2;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site fields (arising from perm. m'poles) //
    // ++++++++++++++++++++++++++++++++++++++++++++++ //

//    cout << endl << "... ... ... 0th-order field" << flush;
    for (sit1 = _polsPolSphere.begin(), seg1 = _segsPolSphere.begin();
         sit1 < _polsPolSphere.end();
         ++sit1, ++seg1) {
    for (sit2 = sit1 + 1, seg2 = seg1 + 1;
         sit2 < _polsPolSphere.end();
         ++sit2, ++seg2) {

        // Intra-pair permanent induction field?
         if ( !induce_intra_pair ) {
             if ( (  (*seg1)->getId() == job->getSeg1Id()
                  || (*seg1)->getId() == job->getSeg2Id() )
               && (  (*seg2)->getId() == job->getSeg1Id()
                  || (*seg2)->getId() == job->getSeg2Id() )) {
                 continue;
             }
         }

         for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
         for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

             _actor.FieldPerm(*(*pit1), *(*pit2));
         }}
    }}

    // +++++++++++++++++++ //
    // 1st-order induction //
    // +++++++++++++++++++ //

//    cout << " | Induce " << endl;
    if (!job->StartFromCPT()) { // OVERRIDE
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

             for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                 (*pit1)->InduceDirect();
             }
        }
    }
    else {
        assert(false); // Load induced dipole moments from file
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

                _actor.FieldIndu(*(*pit1),*(*pit2));                            // <- Intra-pair => zero out
            }}
        }

        // Inter-site contribution to induction field
//        cout << " | Inter-Site (" << iter << ")" << flush;
        //boost::progress_display show_progress( _polsPolSphere.size() );
//        T.restart();

        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//        for (sit1 = _polsPolSphere.begin();
//             sit1 < _polsPolSphere.end();
//             ++sit1) {
//        for (sit2 = sit1 + 1;
//             sit2 < _polsPolSphere.end();
//             ++sit2) {
//
//            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
//            for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {
//
//                _actor.FieldIndu(*(*pit1), *(*pit2));                           // <- Pair-environment => figure sth out
//            }}
//        }}

        for (int id = 0; id < this->_master->_subthreads; ++id) {
            _indus[id]->Start();
        }

        for (int id = 0; id < this->_master->_subthreads; ++id) {
            _indus[id]->WaitDone();
        }

        this->ClearTodoTable();

        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//        cout << " | dt " << T.elapsed() << flush;

        // Induce again
//        cout << " | Induce (" << iter << ")" << flush;
        for (sit1 = _polsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1) {

             for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                 (*pit1)->Induce(wSOR);                                         // <- Okay if alpha = 0
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
            cout << endl << "... ... ... WARNING Induced multipoles for job "
                 << job->getId() << " did not converge to precision: "
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


double XMP::JobXMP::Energy(int state, XJob *job) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    //       PAIR/SITE     <->     SPH1      <->       SPH2        //

    _actor.ResetEnergy();
    double E_Tot = 0.0;
    double E_Pair_Pair = 0.0;   // <- Pair-Pair interaction permanent + induced
    double E_Pair_Sph1 = 0.0;   // <- Pair-Sph1 interaction permanent + induced
    double E_Sph1_Sph1 = 0.0;   // <- Sph1-Sph1 interaction permanent + induced
    double E_Pair_Sph2 = 0.0;   // <- Pair-Sph2 interaction permanent + induced


    vector< Segment* >              ::iterator      seg1;
    vector< Segment* >              ::iterator      seg2;
    vector< vector<PolarSite*> >    ::iterator      sit1;
    vector< vector<PolarSite*> >    ::iterator      sit2;
    vector< PolarSite* >            ::iterator      pit1;
    vector< PolarSite* >            ::iterator      pit2;

    if (job->getType() == "pair") {

        // +++++++++++++++++ //
        // Inter-site energy //
        // +++++++++++++++++ //

        for (sit1 = _polsPolSphere.begin(), seg1 = _segsPolSphere.begin();
             sit1 < _polsPolSphere.end();
             ++sit1, ++seg1) {
        for (sit2 = sit1 + 1, seg2 = seg1 + 1;
             sit2 < _polsPolSphere.end();
             ++sit2, ++seg2) {

            //if ( abs(_top->PbShortestConnect((*seg1)->getPos(),_seg->getPos()))
            //        > _master->_cutoff) { throw runtime_error("Not this."); }

            //cout << "\r... ... Calculating interaction energy for pair "
            //     << (*seg1)->getId() << "|" << (*seg2)->getId() << "   " << flush;

            // Intra-pair interaction?
            if ( ((*seg1)->getId() == job->getSeg1Id() || (*seg1)->getId() == job->getSeg2Id())
              && ((*seg2)->getId() == job->getSeg1Id() || (*seg2)->getId() == job->getSeg2Id()) ) {

                for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

                    //(*pit1)->PrintInfo(cout);
                    //(*pit2)->PrintInfo(cout);

                    E_Pair_Pair += _actor.EnergyInter(*(*pit1), *(*pit2));
                }}
            }

            // Pair-non-pair interaction
            else if ( ((*seg1)->getId() == job->getSeg1Id() || (*seg1)->getId() == job->getSeg2Id())
                    ^ ((*seg2)->getId() == job->getSeg1Id() || (*seg2)->getId() == job->getSeg2Id()) ) {

                for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

                    //(*pit1)->PrintInfo(cout);
                    //(*pit2)->PrintInfo(cout);

                    E_Pair_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }}
            }

            // Non-pair-non-pair interaction?
            else {
                for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

                    //(*pit1)->PrintInfo(cout);
                    //(*pit2)->PrintInfo(cout);

                    E_Sph1_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }}
            }

        }}

        // ++++++++++++++++++ //
        // Outer-Shell energy //
        // ++++++++++++++++++ //

        vector< PolarSite* > central1 = _polarSites_job[ job->getSeg1Id() - 1 ];
        vector< PolarSite* > central2 = _polarSites_job[ job->getSeg2Id() - 1 ];

        for (sit1 = _polsOutSphere.begin(); sit1 < _polsOutSphere.end(); ++sit1) {
            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = central1.begin(); pit2 < central1.end(); ++pit2) {
                    E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));      
                }
            }
        }
        for (sit1 = _polsOutSphere.begin(); sit1 < _polsOutSphere.end(); ++sit1) {
            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = central2.begin(); pit2 < central2.end(); ++pit2) {
                    E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));      
                }
            }
        }


        E_Tot = E_Pair_Pair + E_Pair_Sph1 + E_Sph1_Sph1 + E_Pair_Sph2;

        if (_master->_maverick) {
            cout << endl << "... ... ... ... "
                 << "E(" << state << ") = " << E_Tot * int2eV << " eV "
                 << " = (Pair, pair) " << E_Pair_Pair * int2eV
                 << " + (Pair, Sph1) " << E_Pair_Sph1 * int2eV
                 << " + (Sph1, Sph1) " << E_Sph1_Sph1 * int2eV
                 << " + (Pair, sph2) " << E_Pair_Sph2 * int2eV
                 << flush;
        }


        if (_master->_maverick) {
            cout << endl
                 << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
                 << " eV = (P ~) " << _actor.getEP()    * int2eV
                 << " + (U ~) " << _actor.getEU_INTER() * int2eV
                 << " + (U o) " << _actor.getEU_INTRA() * int2eV
                 << " , with (O ~) " << E_Pair_Sph2 * int2eV << " eV"
                 << flush;
        }

        job->setEnergy(E_Tot*int2eV,           E_Pair_Pair*int2eV,
                       E_Pair_Sph1*int2eV,     E_Sph1_Sph1*int2eV,
                       E_Pair_Sph2*int2eV,
                       _actor.getEP()*int2eV, _actor.getEU_INTER() * int2eV);
    }

    else if (job->getType() == "site") {


        for (int id = 0; id < _master->_subthreads; ++id) {
            _indus[id]->SetSwitch(0);
        }

        // +++++++++++++++++ //
        // Inter-site energy //
        // +++++++++++++++++ //

//        for (sit1 = _polsPolSphere.begin(), seg1 = _segsPolSphere.begin();
//             sit1 < _polsPolSphere.end();
//             ++sit1, ++seg1) {
//        for (sit2 = sit1 + 1, seg2 = seg1 + 1;
//             sit2 < _polsPolSphere.end();
//             ++sit2, ++seg2) {
//
//            // Intra-site interaction?
//            // ... not counted.
//
//            // Site-non-site interaction
//            if ( job->getSiteId() == (*seg1)->getId()
//              || job->getSiteId() == (*seg2)->getId() ) {
//
//                for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
//                for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {
//
//                    //(*pit1)->PrintInfo(cout);
//                    //(*pit2)->PrintInfo(cout);
//
//                    E_Pair_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
//                }}
//            }
//
//            // Non-pair-non-pair interaction?
//            else {
//                for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
//                for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {
//
//                    //(*pit1)->PrintInfo(cout);
//                    //(*pit2)->PrintInfo(cout);
//
//                    E_Sph1_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
//                }}
//            }
//
//        }}

        double eu_inter = 0.0;
        double eu_intra = 0.0;
        double e_perm   = 0.0;

        for (int id = 0; id < this->_master->_subthreads; ++id) {
            _indus[id]->Start();
        }

        for (int id = 0; id < this->_master->_subthreads; ++id) {
            _indus[id]->WaitDone();
        }



        for (int id = 0; id < this->_master->_subthreads; ++id) {
            E_Pair_Pair += _indus[id]->GetEPairPair();
            E_Pair_Sph1 += _indus[id]->GetEPairSph1();
            E_Sph1_Sph1 += _indus[id]->GetESph1Sph1();

            eu_inter += _indus[id]->GetActor().getEU_INTER();
            eu_intra += _indus[id]->GetActor().getEU_INTRA();
            e_perm   += _indus[id]->GetActor().getEP();
        }

        this->ClearTodoTable();






        // ++++++++++++++++++ //
        // Outer-Shell energy //
        // ++++++++++++++++++ //

        vector< PolarSite* > central1 = _polarSites_job[ job->getSeg1Id() - 1 ];

        for (sit1 = _polsOutSphere.begin(); sit1 < _polsOutSphere.end(); ++sit1) {
            for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {
                for (pit2 = central1.begin(); pit2 < central1.end(); ++pit2) {
                    E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));      
                }
            }
        }


        E_Tot = E_Pair_Pair + E_Pair_Sph1 + E_Sph1_Sph1 + E_Pair_Sph2;

        if (_master->_maverick) {
            cout << endl << "... ... ... ... "
                 << "E(" << state << ") = " << E_Tot * int2eV << " eV "
                 << " = (Site, Site) " << E_Pair_Pair * int2eV
                 << " + (Site, Sph1) " << E_Pair_Sph1 * int2eV
                 << " + (Sph1, Sph1) " << E_Sph1_Sph1 * int2eV
                 << " + (Site, Sph2) " << E_Pair_Sph2 * int2eV
                 << flush;
        }


        if (_master->_maverick) {
            cout << endl
                 << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
                 << " eV = (P ~) " << e_perm    * int2eV
                 << " + (U ~) "    << eu_inter  * int2eV
                 << " + (U o) "    << eu_intra  * int2eV
                 << " , with (O ~) " << E_Pair_Sph2 * int2eV << " eV"
                 << flush;
        }

        job->setEnergy(E_Tot*int2eV,           E_Pair_Pair*int2eV,
                       E_Pair_Sph1*int2eV,     E_Sph1_Sph1*int2eV,
                       E_Pair_Sph2*int2eV, 
                       e_perm, eu_inter * int2eV);


        
    }

    else { assert(false); }




    return E_Tot;
}


double XMP::JobXMP::EnergyStatic(int state, XJob *job) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    _actor.ResetEnergy();
    double E_Tot = 0.0;
    double E_Pair_Pair = 0.0;
    double E_Pair_Sph1 = 0.0;
    double E_Sph1_Sph1 = 0.0;
    double E_Pair_Sph2 = 0.0;

    vector< Segment* >              ::iterator      seg1;
    vector< Segment* >              ::iterator      seg2;
    vector< vector<PolarSite*> >    ::iterator      sit1;
    vector< vector<PolarSite*> >    ::iterator      sit2;
    vector< PolarSite* >            ::iterator      pit1;
    vector< PolarSite* >            ::iterator      pit2;

    if (job->getType() == "pair") {

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Interaction pair <-> inner cut-off, without intra-pair interaction //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        vector< PolarSite* > central1 = _polarSites_job[ job->getSeg1Id() - 1 ];
        vector< PolarSite* > central2 = _polarSites_job[ job->getSeg2Id() - 1 ];

        for (seg1 = _segsPolSphere.begin(); seg1 < _segsPolSphere.end(); ++seg1) {

            int id = (*seg1)->getId();

            if (id == job->getSeg1Id() || id == job->getSeg2Id() ) {
                continue;
            }

            for (pit1 = _polarSites_job[id-1].begin();
                 pit1 < _polarSites_job[id-1].end();
                 ++pit1) {
                for (pit2 = central1.begin();
                     pit2 < central1.end();
                     ++pit2) {

                     E_Pair_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
                for (pit2 = central2.begin();
                     pit2 < central2.end();
                     ++pit2) {

                     E_Pair_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
            }
        }


        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Interaction pair <-> outer cut-off                                 //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        for (seg1 = _segsOutSphere.begin(); seg1 < _segsOutSphere.end(); ++seg1) {

            int id = (*seg1)->getId();

            if (id == job->getSeg1Id() || id == job->getSeg2Id() ) {
                throw std::runtime_error("This should not have happened.");
            }

            for (pit1 = _polarSites_job[id-1].begin();
                 pit1 < _polarSites_job[id-1].end();
                 ++pit1) {
                for (pit2 = central1.begin();
                     pit2 < central1.end();
                     ++pit2) {

                     E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
                for (pit2 = central2.begin();
                     pit2 < central2.end();
                     ++pit2) {

                     E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
            }
        }


        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Intra-pair interaction                                             //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        for (pit1 = central1.begin();
             pit1 < central1.end();
             ++pit1) {
        for (pit2 = central2.begin();
             pit2 < central2.end();
             ++pit2) {

            E_Pair_Pair += _actor.EnergyInter(*(*pit1), *(*pit2));
        }}

        E_Tot = E_Pair_Pair + E_Pair_Sph1 + E_Pair_Sph2;

        if (_master->_maverick) {
            cout << endl << "... ... ... ... "
                 << "E(" << state << ") = " << E_Tot * int2eV << " eV "
                 << " = (Pair, intra) " << E_Pair_Pair * int2eV
                 << " + (Pair, inner) " << E_Pair_Sph1 * int2eV
                 << " + (Pair, outer) " << E_Pair_Sph2 * int2eV
                 << flush;
        }

        if (_master->_maverick) {
            cout << endl
                 << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
                 << " eV = (P ~) " << _actor.getEP()       * int2eV
                 << " + (U ~) " << _actor.getEU_INTER() * int2eV
                 << " + (U o) " << _actor.getEU_INTRA() * int2eV << " eV"
                 << ", statics only. "
                 << flush;
        }

        job->setEnergy(E_Tot*int2eV,           E_Pair_Pair*int2eV,
                       E_Pair_Sph1*int2eV,     E_Sph1_Sph1*int2eV,
                       E_Pair_Sph2*int2eV,
                       _actor.getEP()*int2eV, _actor.getEU_INTER() * int2eV);
    }


    else if (job->getType() == "site") {

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Interaction site <-> inner cut-off, without intra-pair interaction //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        vector< PolarSite* > central1 = _polarSites_job[ job->getSiteId() - 1 ];

        for (seg1 = _segsPolSphere.begin(); seg1 < _segsPolSphere.end(); ++seg1) {

            int id = (*seg1)->getId();

            if (id == job->getSiteId()) {
                continue;
            }

            for (pit1 = _polarSites_job[id-1].begin();
                 pit1 < _polarSites_job[id-1].end();
                 ++pit1) {
                for (pit2 = central1.begin();
                     pit2 < central1.end();
                     ++pit2) {

                     E_Pair_Sph1 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
            }
        }


        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Interaction site <-> outer cut-off                                 //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        for (seg1 = _segsOutSphere.begin(); seg1 < _segsOutSphere.end(); ++seg1) {

            int id = (*seg1)->getId();

            if (id == job->getSiteId()) {
                throw std::runtime_error("__ERROR__whx_071");
            }

            for (pit1 = _polarSites_job[id-1].begin();
                 pit1 < _polarSites_job[id-1].end();
                 ++pit1) {
                for (pit2 = central1.begin();
                     pit2 < central1.end();
                     ++pit2) {

                     E_Pair_Sph2 += _actor.EnergyInter(*(*pit1), *(*pit2));
                }
            }
        }


        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
        // Intra-site interaction                                             //
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        // Intra-site energy ...
        // ... not counted.

        E_Tot = E_Pair_Pair + E_Pair_Sph1 + E_Pair_Sph2;

        if (_master->_maverick) {
            cout << endl << "... ... ... ... "
                 << "E(" << state << ") = " << E_Tot * int2eV << " eV "
                 << " = (Site, intra) " << E_Pair_Pair * int2eV
                 << " + (Site, inner) " << E_Pair_Sph1 * int2eV
                 << " + (Site, outer) " << E_Pair_Sph2 * int2eV
                 << flush;
        }

        if (_master->_maverick) {
            cout << endl
                 << "... ... ... ... E(" << state << ") = " << E_Tot * int2eV
                 << " eV = (P ~) " << _actor.getEP()       * int2eV
                 << " + (U ~) " << _actor.getEU_INTER() * int2eV
                 << " + (U o) " << _actor.getEU_INTRA() * int2eV << " eV"
                 << ", statics only. "
                 << flush;
        }

        job->setEnergy(E_Tot*int2eV,           E_Pair_Pair*int2eV,
                       E_Pair_Sph1*int2eV,     E_Sph1_Sph1*int2eV,
                       E_Pair_Sph2*int2eV,
                       _actor.getEP()*int2eV, _actor.getEU_INTER() * int2eV);
    }

    return E_Tot;
}


XMP::XJob *XMP::RequestNextJob(int id, Topology *top) {

    _nextJobMutex.Lock();

    XJob *workOnThis;

    if (_nextXJob == _XJobs.end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextXJob;
        _nextXJob++;
        cout << endl 
             << "... ... Thread " << id << " evaluating job "
             << workOnThis->getId() << " " << workOnThis->getTag()
             << flush;
    }

    _nextJobMutex.Unlock();

    return workOnThis;
}



// ========================================================================== //
//                         XINTERACTOR MEMBER FUNCTIONS                       //
// ========================================================================== //


/**
 * Used in ESP calculator (initialize stage of XMP)
 */
inline double XMP::XInteractor::PotentialPerm(vec r,
                                                     PolarSite &pol) {

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
 * Used in ESF calculator (initialize stage of XMP)
 */
inline vec XMP::XInteractor::FieldPermESF(vec r,
                                                 PolarSite &pol) {
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
inline void XMP::XInteractor::FieldInduAlpha(PolarSite &pol1,
                                                    PolarSite &pol2) {
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _xm->_aDamp;
    u3   = 1 / (R3 * sqrt(pol1.P1 * pol2.P1));

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
inline void XMP::XInteractor::FieldIndu(PolarSite &pol1,
                                               PolarSite &pol2) {

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
    a    = _xm->_aDamp;
    u3   = 1 / (R3 * sqrt(pol1.P1 * pol2.P1));

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
inline void XMP::XInteractor::FieldPerm(PolarSite &pol1,
                                               PolarSite &pol2) {

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
inline double XMP::XInteractor::EnergyIntra(PolarSite &pol1,
                                                   PolarSite &pol2) {    

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;
    
    // Thole damping init.
    a    = _xm->_aDamp;
    u3   = 1 / (R3 * sqrt(pol1.P1 * pol2.P1));

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
inline double XMP::XInteractor::EnergyInter(PolarSite &pol1,
                                                   PolarSite &pol2) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _xm->_aDamp;
    u3   = 1 / (R3 * sqrt(pol1.P1 * pol2.P1));


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
inline double XMP::XInteractor::EnergyInterESP(PolarSite &pol1,
                                                      PolarSite &pol2) {
    
    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = pol2.getPos() - pol1.getPos();
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

    // Thole damping init.
    a    = _xm->_aDamp;
    u3   = 1 / (R3 * sqrt(pol1.P1 * pol2.P1));

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



}}

#endif

