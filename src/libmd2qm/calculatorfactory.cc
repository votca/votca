#include "calculatorfactory.h"

#include "calculators/integrals.h"
#include "calculators/marcusrates.h"
#include "calculators/readxml.h"
#include "calculators/writexml.h"
#include "calculators/contkmc.h"
#include "calculators/estatics.h"


void CalculatorFactory::RegisterAll(void)
{
        Calculators().Register<CalcIntegrals>("integrals");
        Calculators().Register<WriteXML>("writexml");
        Calculators().Register<ReadXML>("readxml");
        Calculators().Register<ContKmc>("kmc");
        Calculators().Register<CalcEstatics>("estat");
        Calculators().Register<MarcusRates>("marcusrates");
        Calculators().Register<CalcIntegrals>("histintegrals");
}