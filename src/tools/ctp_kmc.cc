/* 
 * File:   main.cpp
 * Author: vehoff
 *
 * Created on April 21, 2010, 10:32 AM
 */

#include <stdlib.h>

#include "qmapplication.h"
#include "cont_kmc.h"

/*
 * 
 */
int main(int argc, char** argv) {

    QMApplication qm_app;
    ContKmc cont_kmc;

    qm_app.AddCalculator(dynamic_cast<QMCalculator*>(&cont_kmc));

    qm_app.Run(argc, argv);

    return (EXIT_SUCCESS);
}

