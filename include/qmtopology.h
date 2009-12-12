/* 
 * File:   qmtopology.h
 * Author: james
 *
 * Created on 12 December 2009, 13:14
 */

#ifndef _QMTOPOLOGY_H
#define	_QMTOPOLOGY_H

#include <votca/csg/topology.h>
#include <votca/csg/nblist.h>

/**
    \brief contains the atomistic representation of the QM beads

 * It contains the QMBeads.
*/

class QMTopology : Topology{
public:
    QMTopology();
    ~QMTopology();
    int CreateBeads();
    int Init();
private:
    CreateQMBeads();

    NBList _neighs;
    map <> ;
};

#endif	/* _QMTOPOLOGY_H */

