#include <votca/ctp/calculatorfactory.h>
#include "votca_config.h"

#include "calculators/sandbox.h"
#include "calculators/neighborlist.h"
#include "calculators/stateserver.h"
#include "calculators/tdump.h"
#include "calculators/rates.h"
#include "calculators/izindo.h"
#include "calculators/einternal.h"
#include "calculators/eoutersphere.h"
#include "calculators/emultipole.h"
#include "calculators/emultipole_stdal.h"




namespace votca { namespace ctp {

void Calculatorfactory::RegisterAll(void)
{	
        Calculators().Register<Sandbox>("sandbox"); // Test calculator
        Calculators().Register<Neighborlist>("neighborlist");
        Calculators().Register<StateServer>("stateserver");
        Calculators().Register<TDump>("tdump");
        Calculators().Register<IZindo>("izindo");
        Calculators().Register<Rates>("rates");
        Calculators().Register<EInternal>("einternal");
        Calculators().Register<EOutersphere>("eoutersphere");
        Calculators().Register<EMultipole>("emultipole");
        Calculators().Register<EMultipole_StdAl>("emultipole2");
}

}}
