/* 
 * File:   EMultipole.h
 * Author: poelking
 *
 * Created on December 27, 2011, 4:10 PM
 */

#ifndef EMULTIPOLE_H
#define	EMULTIPOLE_H

#include <votca/ctp/qmcalculator2.h>
#include <votca/ctp/topology.h>

namespace votca { namespace ctp {


class EMultipole : public QMCalculator2
{
public:

    EMultipole() { };
   ~EMultipole() { };

    void Initialize(Topology *top, Property *opt);
    void GetMPoles(string &xmlfile);
    void EquipTop(Topology *top);

    bool EvaluateFrame(Topology *top);
    void ChargeMol(Molecule *mol, int state);
    void Induce(Topology *top);
    void Depolarize(Topology *top);

    void CalcIntEnergy(Topology *top, int &state);
    double PairIntEnergy(Atom *atm, Atom *btm);

    string Identify() { return "Thole Calculator"; }
    void PrintInfo(const string &key);
    void PrintProgBar(int percent);




private:

    // interaction parameters
    bool            _useCutoff;
    double          _cutoff;
    bool            _useExp;
    double          _aDamp;
    bool            _useScaling;
    vector<double>  _scale1;

    // convergence parameters
    float           _omegSOR;
    double          _epsTol;
    int             _maxIter;

    // store multipole data (polarizability + charges);
    // keys are atom identifiers "atmname_rsdname_molname"
    // TODO Move to atom.h or similar to get rid
    // TODO of map structure and increase look-up speed
    map < string, matrix >           _ptensors;
    map < string, matrix >::iterator _ptit;
    map < string, map<int, double> >           _chrgs;
    map < string, map<int, double> >::iterator _chrgit;
    map < string, map<int, bool > >            _hasChrg;

    // occupations for which to calculate site energies;
    // neutral assumed as precondition
    bool _anion;
    bool _cation;


    // store induced dipoles for atoms in topology;
    // order is the same as in top-> atom container
    vector <vec> _muInd;
    // field at sites due to permanent multipoles
    vector <vec> _pField;
    // field at sites due to permanent multipoles
    vector <vec> _iField;
    // TODO Move polarization tensor to Atom class
    vector <matrix> _dpolT;

    // store site-energy contributions;
    // keys are -1 ("anion"), 0 ("neutral"), 1 ("cation");
    // re-defined for each segment in ::CalcIntEnergy;
    map < int, double > _pInter; // permanent
    map < int, double > _iInter; // induction
    map < int, double > _pIntra; // permanent
    map < int, double > _iIntra; // induction

    // Thole damping functions
    inline double lambda3(double &u) {
        return 1 - exp( -1*_aDamp*pow(u,3));
    }
    inline double lambda5(double &u) {
        return 1 - ( 1 + _aDamp*pow(u,3) ) *
                     exp( -1*_aDamp*pow(u,3) );
    }

    // some interaction tensors
    inline double T0(vec &R) { return 1 / abs(R); }

    inline vec T1(vec &R) { return -1/pow(abs(R),3) * R; }
    inline vec T1(vec &R, double & u) { return -1/pow(abs(R),3)*lambda3(u)*R; }

    inline matrix T2(vec &R) {
        matrix diag = matrix();
        diag.UnitMatrix();
        diag *= 1/pow(abs(R),3);
        return ( 3/pow(abs(R),5) * R|R ) - diag;
    }
    inline matrix T2(vec &R, double u) {
        matrix diag = matrix();
        diag.UnitMatrix();
        diag *= 1/pow(abs(R),3) * lambda3(u);
        return ( 3/pow(abs(R),5) * lambda5(u) * R|R ) - diag;
    }

};




}} /* exit namespaces votca, ctp */

#endif /* EMULTIPOLE_H */