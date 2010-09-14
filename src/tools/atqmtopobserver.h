/* 
 * File:   atqmtopobserver.h
 * Author: james
 *
 * Created on March 10, 2010, 9:52 AM
 */

#ifndef _ATQMTOPOBSERVER_H
#define	_ATQMTOPOBSERVER_H

#include "qmapplication.h"
#include <votca/csg/trajectorywriter.h>

/// will read in the atomistic trajectory and
/// generate corresponding coarse graines
/// and qm topologies according to the list charge definition
class AtQmObserver: public QMApplication
{
public:
    AtQmObserver();
    ~AtQmObserver();

    string ProgramName() { return "ctp_printqmattop"; }
    void HelpText(ostream &out) {
        out << "dump results of mapping procedure";
    }

    void Initialize();
    bool EvaluateOptions();

    bool EvaluateFrame();
    void EndEvaluate();

protected:
    TrajectoryWriter *_writerCG;
    TrajectoryWriter *_writerQM;
};

#endif	/* _ATQMTOPOBSERVER_H */

