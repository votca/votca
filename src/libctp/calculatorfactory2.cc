#include <votca/ctp/calculatorfactory2.h>
#include "votca_config.h"

#include "calculators/sandbox2.h"
#include "calculators/neighborlist2.h"
#include "calculators/stateserver.h"
#include "calculators/tdump2.h"
#include "calculators/rates2.h"
#include "calculators/izindo.h"
#include "calculators/emultipole.h"
#include "calculators/emultipole2.h"
#include "calculators/emultipole3.h"



namespace votca { namespace ctp {

void CalculatorFactory2::RegisterAll(void)
{	
        Calculators().Register<Sandbox2>("sandbox"); // Test calculator
        Calculators().Register<Neighborlist2>("neighborlist");
        Calculators().Register<StateServer>("stateserver");
        Calculators().Register<TDump>("tdump");
        Calculators().Register<IZindo>("izindo");
        Calculators().Register<Rates2>("rates");
        Calculators().Register<EMultipole>("emultipole_single");
        Calculators().Register<EMultipole2>("emultipole");
        Calculators().Register<EMultipole3>("emultipole3");
}

}}
