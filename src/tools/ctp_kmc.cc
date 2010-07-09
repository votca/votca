/* 
 * File:   main.cpp
 * Author: vehoff
 *
 * Created on April 21, 2010, 10:32 AM
 */

#include <stdlib.h>

#include "qmapplication.h"
#include "calculatorfactory.h"

/*
 * 
 */
int main(int argc, char** argv) {

    QMApplication qm_app;

    qm_app.AddCalculator(Calculators().Create("kmc"));

    qm_app.Run(argc, argv);

    return (EXIT_SUCCESS);
}

