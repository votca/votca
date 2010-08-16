#include "calculatorfactory.h"

#include "calculators/integrals.h"
#include "calculators/marcusrates.h"
#include "calculators/readxml.h"
#include "calculators/writexml.h"
#include "calculators/contkmc.h"
#include "calculators/estatics.h"
#include "calculators/histintegrals.h"
#include "calculators/shufflenrg.h"
#include "calculators/generate_nrgs.h"
#include "calculators/energycorr.h"
#include "calculators/polymerrates.h"

void CalculatorFactory::RegisterAll(void)
{
        Calculators().Register<CalcIntegrals>("integrals");
        Calculators().Register<WriteXML>("writexml");
        Calculators().Register<ReadXML>("readxml");
        Calculators().Register<ContKmc>("kmc");
        Calculators().Register<CalcEstatics>("estat");
        Calculators().Register<MarcusRates>("marcusrates");
        Calculators().Register<CalcHistIntegrals>("histintegrals");
        Calculators().Register<ShuffleNrg>("shufflenrg");
        Calculators().Register<GenerateNrgs>("generatenrgs");
        Calculators().Register<EnergyCorr>("energycorr");
        Calculators().Register<PolymerRates>("polymerrates");
}