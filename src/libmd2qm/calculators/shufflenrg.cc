/* 
 * File:   shufflenrg.cc
 * Author: lukyanov
 * 
 * Created on July 21, 2010, 2:48 PM
 */

#include <stdlib.h>
#include "shufflenrg.h"
#include <math.h>
#include <list>

bool ShuffleNrg::EvaluateFrame(QMTopology *top) {

    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;
    vector<double> energies;
    int i_en = 0;

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        energies.push_back( (*itl)->getEnergy() );
    }
    
    random_shuffle ( energies.begin(), energies.end() );

    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
        (*itl)->setEnergy( energies[i_en] );
        i_en++;
    }


    return true;
}