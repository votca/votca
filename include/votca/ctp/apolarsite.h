#ifndef APOLARSITE_H
#define APOLARSITE_H

#include <cstdlib>
#include <fstream>
#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>


namespace votca { namespace ctp {


using namespace votca::tools;


class Topology;
class Molecule;
class Segment;
class Fragment;
class QMThread;


class APolarSite
{

    friend class BasicInteractor;
    friend class MolPol;
    friend class ZMultipole;
    friend class QMultipole;
    friend class XInteractor;
    friend class QMMIter;


public:

    APolarSite(int id, string name)
            : _id(id), _name(name), _isVirtual(false), _locX(vec(1,0,0)),
              _locY(vec(0,1,0)),    _locZ(vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

    APolarSite()
            : _id(-1),  _isVirtual(false), _locX(vec(1,0,0)),
              _locY(vec(0,1,0)), _locZ(vec(0,0,1))
            { _Qs.resize(3); _Ps.resize(3); };

   ~APolarSite() {};

    int             &getId() { return _id; }
    string          &getName() { return _name; }
    vec             &getPos() { return _pos; }
    int             &getRank() { return _rank; }
    Topology        *getTopology() { return _top; }
    Segment         *getSegment() { return _seg; }
    Fragment        *getFragment() { return _frag; }


    void            setIsVirtual(bool isVirtual) { _isVirtual = isVirtual; }
    void            setPos(vec &pos) { _pos = pos; }
    void            setRank(int rank) { _rank = rank; } // rank; } // OVERRIDE
    void            setTopology(Topology *top) { _top = top; }
    void            setSegment(Segment *seg) { _seg = seg; }
    void            setFragment(Fragment *frag) { _frag = frag; }
    
    bool            getIsVirtual() { return _isVirtual; }
    vector<double> &getQs(int state) { return _Qs[state+1]; }
    void            setQs(vector<double> Qs, int state) { _Qs[state+1] = Qs; }
    void            setPs(matrix polar, int state) { _Ps[state+1] = polar; }
    matrix         &getPs(int state) { return _Ps[state+1]; }
    double          getIsoP() { return 1./3. * (Pxx+Pyy+Pzz); }
    double          getProjP(vec &dir);
    
    double          setQ00(double q, int s) { Q00 = q; _Qs[s+1][0] = q; }
    double         &getQ00() { return Q00; }
    void            Charge(int state);
    void            ChargeDelta(int state1, int state2);

    void            Induce(double wSOR = 0.25);
    void            InduceDirect();
    void            ResetFieldU() { FUx = FUy = FUz = 0.0; }
    void            ResetFieldP() { FPx = FPy = FPz = 0.0; }
    void            ResetU1Hist() { U1_Hist.clear(); }
    void            Depolarize();
    double          HistdU();


    void            ImportFrom(APolarSite *templ, string tag = "basic");
    void            Translate(const vec &shift);
    void            Rotate(const matrix &rot, const vec &refPos);

    void            PrintInfo(std::ostream &out);
    void            PrintTensorPDB(FILE *out, int state);
    void            WriteChkLine(FILE *, vec &, bool, string, double);
    void            WriteXyzLine(FILE *, vec &, string);
    void            WritePdbLine(FILE *out, const string &tag = "");

    vector<APolarSite*> CreateFrom_MPS(string filename, int state) { ; }




private:

    int     _id;
    string  _name;
    bool    _isVirtual;
    vec     _pos;
    vec     _locX;
    vec     _locY;
    vec     _locZ;

    Topology *_top;
    Segment  *_seg;
    Fragment *_frag;

    vector < vector<double> > _Qs;
    int     _rank;

    vector < matrix > _Ps;
    
    double Pxx, Pxy, Pxz;                   // Dipole polarizability tensor
    double      Pyy, Pyz;
    double           Pzz;

    vec pax; double eigenpxx;               // Principal axes, eigenvals. of P
    vec pay; double eigenpyy;
    vec paz; double eigenpzz;

    double eigendamp;

    double Q00;
    double Q1x, Q1y, Q1z;
    double Q20, Q21c, Q21s, Q22c, Q22s;

    double U1x, U1y, U1z;                   // Induced dipole
    double FPx, FPy, FPz;                   // Electric field (due to permanent)
    double FUx, FUy, FUz;                   // Electric field (due to induced)
    vector< vec > U1_Hist;                  // Ind. u history





};




vector<APolarSite*> APS_FROM_MPS(string filename, int state, QMThread *thread = NULL);





class BasicInteractor
{
public:

    BasicInteractor() {};
   ~BasicInteractor() {};

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

    void        SetADamp(double aDamp) { a = aDamp; }

    inline void FieldInduAlpha(APolarSite &pol1, APolarSite &pol2)  {
        // NOTE >>> Exp. damping factor a must have been set     <<< NOTE //
        // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
        e12  = pol2.getPos() - pol1.getPos();
        R    = 1/abs(e12);
        R2   = R*R;
        R3   = R2*R;
        R4   = R3*R;
        R5   = R4*R;
        e12 *= R;

        // Thole damping init.
//        double P_proj_1 = fabs(e12 * pol1.pax) * pol1.eigenpxx
//                        + fabs(e12 * pol1.pay) * pol1.eigenpyy
//                        + fabs(e12 * pol1.paz) * pol1.eigenpzz;
//        double P_proj_2 = fabs(e12 * pol2.pax) * pol2.eigenpxx
//                        + fabs(e12 * pol2.pay) * pol2.eigenpyy
//                        + fabs(e12 * pol2.paz) * pol2.eigenpzz;

        //u3   = 1 / (R3 * sqrt(pol1.getIsoP() * pol2.getIsoP()));
        u3   = 1 / (R3 * sqrt(pol1.getProjP(e12) * pol2.getProjP(e12)));

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


private:

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

};



}}


#endif
