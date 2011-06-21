#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <stdexcept>
#include "atqmtopobserver.h"

int main(int argc, char** argv)
{
    AtQmObserver int_calc;
    return int_calc.Exec(argc, argv);
}
