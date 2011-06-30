#ifndef _DUMP_TRAJECORY_H
#define	_DUMP_TRAJECTORY_H

#include <stdlib.h>
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

#include "qmpair.h"
#include "qmcalculator.h"

#include <votca/csg/trajectorywriter.h>

/**
    \brief Outputs the coarse-grained and back-mapped (using rigid fragments) trajectories


Callname: tdump

Useful for checking whether the mapping of the atomistic trajectory on conjugated segments and rigid fragments is correct. One can use VisualMolecularDnamics (vmd) to view the initial, coarse-grained, and back-mapped trajectories together.

*/
class DumpTrajectory : public QMCalculator
{
public:
    DumpTrajectory() {};
    ~DumpTrajectory() {};

    const char *Description() { return "Outputs the coarse-grained and back-mapped (using rigid fragments) trajectories"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);
   
private:
    Property * _options;
    string _nameCG, _nameQM;
    TrajectoryWriter *_writerCG, *_writerQM; 
};

inline void DumpTrajectory::Initialize(QMTopology *top, Property *options) {
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

inline bool DumpTrajectory::EvaluateFrame(QMTopology *top) {
    
     // dumping the coarse-grained trajectory
    _writerCG->Write( top );
    
    // creating the back-mapped atomistic trajectory
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    Topology qmAtomisticTop;
    for (vector<QMCrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++){
        top->AddAtomisticBeads(*itl,&qmAtomisticTop);
    }
    
    // dumping the back-mapped trajectory and cleaning
    _writerQM->Write(&qmAtomisticTop);
    
    qmAtomisticTop.Cleanup();

    return true;
}

inline void DumpTrajectory::EndEvaluate(QMTopology *top)
{
    _writerCG->Close();
    _writerQM->Close();
}

#endif	/* _DUMP_TRAJECTORY_H */

