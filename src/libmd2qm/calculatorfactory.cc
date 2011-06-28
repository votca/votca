#include "calculatorfactory.h"

#include "votca_config.h"
#include "calculators/integrals.h"
#include "calculators/marcusrates.h"
#include "calculators/readxml.h"
#include "calculators/writexml.h"
#include "calculators/ecoulomb.h"
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
#include "calculators/sqlitewriter.h"
#include "calculators/nbgen.h"
#include "calculators/occequilibrium.h"
#include "calculators/avgvelocity.h"
#include "calculators/lambdaout.h"
#include "calculators/dump_atoms.h"
#include "calculators/dump_trajectory.h"
#include "calculators/dump_atoms_bj.h"

#ifdef WITH_VOTCA_KMCOLD        
    #include "calculators/contkmc.h"
#endif

void CalculatorFactory::RegisterAll(void)
{
Calculators().Register<CalcIntegrals>("integrals");
        Calculators().Register<WriteXML>("writexml");
        Calculators().Register<ReadXML>("readxml");
        Calculators().Register<CalcEstatics>("ecoulomb");
        Calculators().Register<CalcLambdaOut>("lambdaout");
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
        Calculators().Register<SQLiteWriter>("sqlitewriter");
        Calculators().Register<NBGen>("nbgen");
        Calculators().Register<OccEquilibrium>("occequilibrium");
        Calculators().Register<AvgVelocity>("avgvelocity");
        Calculators().Register<DumpAtomsBJ>("dumpatomsbj");
        Calculators().Register<DumpTrajectory>("dumptraj");
        Calculators().Register<DumpTrajectory>("dumpatomsbj");
#ifdef WITH_VOTCA_KMCOLD        
        Calculators().Register<ContKmc>("kmc");
#endif

}
