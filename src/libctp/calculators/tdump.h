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

#ifndef _TDUMP_H
#define	_TDUMP_H

#include <stdlib.h>
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmcalculator.h>

#include <votca/csg/trajectorywriter.h>

namespace votca { namespace ctp {

/**
    \brief Outputs the coarse-grained and back-mapped (using rigid fragments) trajectories


Callname: tdump

Useful for checking whether the mapping of the atomistic trajectory on conjugated segments and rigid fragments is correct. One can use VisualMolecularDnamics (vmd) to view the initial, coarse-grained, and back-mapped trajectories together.

*/
class Tdump : public QMCalculator
{
public:
    Tdump() {};
    ~Tdump() {};

    const char *Description() { return "Outputs the coarse-grained and back-mapped (using rigid fragments) trajectories"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);
   
private:
    Property * _options;
    string _nameCG, _nameQM;
    TrajectoryWriter *_writerCG, *_writerQM; 
};

inline void Tdump::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if ( options->exists("options.tdump.cg") && 
          options->exists("options.tdump.qm") ) {
        _nameCG = options->get("options.tdump.cg").as<string>();
        _nameQM = options->get("options.tdump.qm").as<string>();
        cout << "Writing the  conjugated  segments trajectory to " << _nameCG <<endl;
        cout << "Writing the backmapped atomistic trajectory to " << _nameQM <<endl;
    } else {
        _nameCG="traj_cg.pdb";
        _nameQM="traj_qm.pdb";
        cout << "Warning: at least one of the trajectory names was not given" << endl
                << "Using default names traj_cg.pdb and traj_qm.pdb" << endl;
    }
    
    string extCG  = _nameCG.substr(_nameCG.length()-4,4);
    string extQM  = _nameCG.substr(_nameQM.length()-4,4);

    _writerCG = TrjWriterFactory().Create(_nameCG);
    if(_writerCG == NULL) throw runtime_error(string("output format not supported: ")+ extCG);

    _writerQM = TrjWriterFactory().Create(_nameQM);
    if(_writerQM == NULL) throw runtime_error(string("output format not supported: ")+ extQM);

    _writerCG->Open(_nameCG);
    _writerQM->Open(_nameQM);

}

inline bool Tdump::EvaluateFrame(QMTopology *top) {
    
     // dumping the coarse-grained trajectory
    _writerCG->Write( top );
    
    // creating the back-mapped atomistic trajectory
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    Topology qmAtomisticTop;
    qmAtomisticTop.setBox(top->getBox());
    qmAtomisticTop.setStep(top->getStep());
    qmAtomisticTop.setTime(top->getTime());


    for (vector<QMCrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++){
        top->AddAtomisticBeads(*itl,&qmAtomisticTop);
    }
    
    // dumping the back-mapped trajectory and cleaning
    _writerQM->Write(&qmAtomisticTop);
    
    qmAtomisticTop.Cleanup();

    return true;
}

inline void Tdump::EndEvaluate(QMTopology *top)
{
    _writerCG->Close();
    _writerQM->Close();
}

}}

#endif	/* _TDUMP_H */

