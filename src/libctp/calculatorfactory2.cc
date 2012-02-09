#include <votca/ctp/calculatorfactory2.h>
#include "votca_config.h"

#include "calculators/sandbox2.h"
#include "calculators/neighborlist2.h"
#include "calculators/stateserver.h"
#include "calculators/tdump2.h"
#include "calculators/emultipole.h"


namespace votca { namespace ctp {

void CalculatorFactory2::RegisterAll(void)
{	
        Calculators().Register<Sandbox2>("sandbox"); // Test calculator
        Calculators().Register<Neighborlist2>("neighborlist");
        Calculators().Register<StateServer>("stateserver");
        Calculators().Register<TDump>("tdump");
        Calculators().Register<EMultipole>("emultipole");
}

}}
