/* 
 * File:   polymer.h
 * Author: james
 *
 * Created on March 29, 2010, 1:24 PM
 */

#ifndef _POLYMER_H
#define	_POLYMER_H

#include "qmapplication.h"
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>



class Polymer : public QMApplication
{
public:
    Polymer();
    ~Polymer();

    void HelpText();
    void AddSpecificOptions();
    void Initialize();
    bool EvaluateFrame();
    void EndEvaluate();

private:
    
    struct WaveFunction{
        double _nrg;
        gsl_vector * _wf;
        int _molid;
        int _wfid;
        int _id;
        vec _pos;
    };
    typedef pair <WaveFunction *, WaveFunction *> PairWF ;
    void UpdatePolTop();
    void UpdateJs(CrgUnit*, CrgUnit *, const double & );
    void UpdatedR(const PairWF & );
    //double ComputeDj(Molecule * one, Molecule *two, WaveFunction *a, WaveFunction *b,const double & J);
    void CalcWaveFunction(Molecule * mol);
    
    void Clean();
    void Save(string &);
    
    /// energy of each wavefunction and transfer integrals
    double _cutnrg;
    vector < WaveFunction * > _states;
    map < PairWF , double > _polJs;
    map < PairWF, vec > _poldR;
    map < int, int> _mcrg2bs;
    vector < vector < WaveFunction * > *> _mstates;

};
/*
inline double dot (gsl_vector * a, gsl_vector * b){
    double r=0.;
    for (int i=0;i< a->size;++i){
        r += gsl_vector_get(a,i) * gsl_vector_get(b,i);
    }
    return r;
}*/
#endif	/* _POLYMER_H */

