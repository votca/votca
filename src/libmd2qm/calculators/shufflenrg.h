#ifndef _SHUFFLENRG_H
#define	_SHUFFLENRG_H

#include "qmpair.h"
#include "qmcalculator.h"

/** \brief Shuffle site energies to get rid of spatial correlations, while preserving the distribution width 

The site energies are randomly assigned to sites which keeps the distribution of site energies but destroys the spatial correlations present e.g. due to charge-dipole interactions.
Therefore the new distribution of site-energy differences may differ from the old distribution of site-energy differences.

Callname: 
*/
class ShuffleNrg : public QMCalculator
{
public:
    ShuffleNrg() {};
    ~ShuffleNrg() {};


    bool EvaluateFrame(QMTopology *top);

};

#endif	/* _SHUFFLENRG_H */

