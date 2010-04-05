#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <stdexcept>
#include "atqmtopobserver.h"

using namespace std;


int main(int argc, char** argv)
{
    AtQmObserver int_calc;
    int_calc.Run(argc, argv);

    return (EXIT_SUCCESS);
}



