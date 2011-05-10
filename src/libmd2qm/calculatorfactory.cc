#include "calculatorfactory.h"

#include "calculators/integrals.h"
#include "calculators/marcusrates.h"
#include "calculators/readxml.h"
#include "calculators/writexml.h"
#include "calculators/contkmc.h"
#include "calculators/estatics.h"
#include "calculators/histintegrals.h"
#include "calculators/histenergeticdisorder.h"
#include "calculators/shufflenrg.h"
#include "calculators/generate_nrgs.h"
#include "calculators/energycorr.h"
#include "calculators/polymerrates.h"
#include "calculators/pairdump.h"
#include "calculators/CurrentDensity.h"
#include "calculators/jortnerrates.h"
#include "calculators/marcusrateslambdaouter.h"
#include "calculators/hdf5writer.h"
#include "calculators/sqlitewriter.h"
#include "calculators/nbgen.h"
#include "calculators/occequilibrium.h"
#include "calculators/avgvelocity.h"
 
void CalculatorFactory::RegisterAll(void)
{
        Calculators().Register<CalcIntegrals>("integrals");
        Calculators().Register<WriteXML>("writexml");
        Calculators().Register<ReadXML>("readxml");
        Calculators().Register<ContKmc>("kmc");
        Calculators().Register<CalcEstatics>("estat");
        Calculators().Register<MarcusRates>("marcusrates");
        Calculators().Register<CalcHistIntegrals>("histintegrals");
        Calculators().Register<CalcHistEnergeticDisorder>("histenergeticdisorder");
        Calculators().Register<ShuffleNrg>("shufflenrg");
        Calculators().Register<GenerateNrgs>("generatenrgs");
        Calculators().Register<EnergyCorr>("energycorr");
        Calculators().Register<PolymerRates>("polymerrates");
        Calculators().Register<PairDump>("pairdump");
        Calculators().Register<CurrentDensity>("currden");
        Calculators().Register<JortnerRates>("jortnerrates");
        Calculators().Register<MarcusRatesLambdaOuter>("marcusrateslambdaouter");
        Calculators().Register<HDF5Writer>("hdf5writer");
        Calculators().Register<SQLiteWriter>("sqlitewriter");
        Calculators().Register<NBGen>("nbgen");
        Calculators().Register<OccEquilibrium>("occequilibrium");
        Calculators().Register<AvgVelocity>("avgvelocity");
}
