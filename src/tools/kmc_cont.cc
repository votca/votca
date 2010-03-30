/* 
 * File:   main.cpp
 * Author: vehoff
 *
 * Created on March 4, 2010, 4:43 PM
 */

#include <stdlib.h>
#include "kmc_cont_app.h"

/*
 * 
 */
int main(int argc, char** argv) {

    KmcCont kmc_cont;
    kmc_cont.Run(argc, argv);

    return (EXIT_SUCCESS);
}

