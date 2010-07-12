/* 
 * File:   ctp_run.cc
 * Author: schrader
 *
 * Created on July 12, 2010, 3:11 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include <calculatorfactory.h>
/*
 *
 */
int main(int argc, char** argv) {

//<<<<<<< local
    QMApplication qm_app;


    qm_app.AddCalculator(Calculators().Create("integrals"));
    qm_app.AddCalculator(Calculators().Create("histintegrals"));


    //    qm_app.AddCalculator(Calculators().Create("writexml"));

    qm_app.Run(argc, argv);

//=======
//    QMApplication calc_estatics;
//    CalcEstatics calc_estat;
//    calc_estatics.AddCalculator(dynamic_cast<QMCalculator*>(&calc_estat));
//    calc_estatics.Run(argc, argv);
//
//>>>>>>> other
    return (EXIT_SUCCESS);
}

