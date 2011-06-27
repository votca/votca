#ifndef _DUMP_TRAJECORY_H
#define	_DUMP_TRAJECTORY_H

#include "qmpair.h"
#include "qmcalculator.h"
#include <votca/csg/trajectorywriter.h>

/**
    \brief Outputs the coarse-grained and back-mapped (using rigid fragments) trajectories


Callname: dumptraj

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

#endif	/* _DUMP_TRAJECTORY_H */

