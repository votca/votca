/* 
 * File:   main.cc
 * Author: vehoff
 *
 * Created on February 2, 2010, 11:25 AM
 */

#include <stdlib.h>
#include "ratecalcapp.h"

/*
 * 
 */
int main(int argc, char** argv) {

    RateCalculator rate_calc;
    rate_calc.Run(argc, argv);

    return (EXIT_SUCCESS);
}

