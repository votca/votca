/* 
 * File:   QMTopologyCreator.h
 * Author: ruehle
 *
 * Created on July 12, 2011, 4:28 PM
 */

#ifndef __VOTCA_CTP_QMTOPOLOGYCREATOR_H
#define	__VOTCA_CTP_QMTOPOLOGYCREATOR_H

#include "qmtopology.h"
class QMTopologyCreator
: public QMTopology
{
public:
    /// update the topology based on cg positons
    void Initialize(Topology &cg_top);
    ///Initialises the charge units
    void InitChargeUnits();

};


#endif	/* QMTOPOLOGYCREATOR_H */

