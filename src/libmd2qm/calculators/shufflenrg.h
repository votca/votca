/* 
 * File:   shufflenrg.h
 * Author: lukyanov
 *
 * Created on July 21, 2010, 2:48 PM
 */

#ifndef _SHUFFLENRG_H
#define	_SHUFFLENRG_H

#include "qmpair.h"
#include "qmcalculator.h"

/// Shuffle array of site energies to get rid of spatial correlations, while preserving the distribution

class ShuffleNrg : public QMCalculator
{
public:
    ShuffleNrg() {};
    ~ShuffleNrg() {};


    bool EvaluateFrame(QMTopology *top);

};

#endif	/* _SHUFFLENRG_H */

