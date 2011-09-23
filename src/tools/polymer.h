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

#ifndef _POLYMER_H
#define	_POLYMER_H

#include <votca/ctp/qmapplication.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace votca::ctp;


class Polymer : public QMApplication
{
public:
 
    void HelpText() {}

    string ProgramName() { return "ctp_polymer"; }
    void HelpText(std::ostream &out) {
        out << "Charge transport in conjugated polymers [OUTDATED]" << endl ;
    }
    
    void Initialize();
    bool EvaluateOptions();

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

