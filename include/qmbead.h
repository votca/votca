/* 
 * File:   QMBead.h
 * Author: james
 *
 * Created on 12 December 2009, 13:23
 */

#ifndef _QMBEAD_H
#define	_QMBEAD_H

#include "moo/crgunit.h"

/**
    \brief contains the position, etc. information for the atomistic representation of 
 * a crgunit

 * The Bead class describes an crgunit. It belongs to a QMTopology.
 * Since it is overloaded from Bead it can use the votca neighbour list codes etc
*/


class QMBead:Bead{
public:
    QMBead(BeadContainer &, CrgUnit &);
    ~QMBead();
    
private:

    /// the energy of the site
    double _nrg;
    
};

#endif	/* _QMBEAD_H */

