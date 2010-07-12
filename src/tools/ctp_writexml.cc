/* 
 * File:   playground.cpp
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:21 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include "calculatorfactory.h"

/*
 *
 */
int main(int argc, char** argv) {

    QMApplication qm_app;

    qm_app.AddCalculator(Calculators().Create("writexml"));
   
    qm_app.Run(argc, argv);

    return (EXIT_SUCCESS);
}

