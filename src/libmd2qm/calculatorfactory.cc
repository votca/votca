#include <votca/ctp/calculatorfactory.h>

#include "votca_config.h"
#include "calculators/ecoulomb.h"
#include "calculators/ehistogram.h"
#include "calculators/egaussian.h"
#include "calculators/ecorrelation.h"
#include "calculators/eshuffle.h"

#include "calculators/izindo.h"
#include "calculators/ihistogram.h"

#include "calculators/tdump.h"

#include "calculators/vaverage.h"

#include "calculators/oboltzmann.h"

#include "calculators/neighborlist.h"
#include "calculators/etinker.h"


#include "calculators/marcusrates.h"
//#include "calculators/readxml.h"
#include "calculators/writexml.h"
#include "calculators/polymerrates.h"
#include "calculators/pairdump.h"
#include "calculators/CurrentDensity.h"
#include "calculators/jortnerrates.h"
#include "calculators/marcusrateslambdaouter.h"
#include "calculators/sqlitewriter.h"
#include "calculators/lambdaout.h"
#include "calculators/dump_atoms_bj.h"
#include "calculators/rates.h"

#ifdef WITH_VOTCA_KMCOLD        
    #include "calculators/contkmc.h"
#endif

void CalculatorFactory::RegisterAll(void)
{
	Calculators().Register<CalcIntegrals>("izindo");
        Calculators().Register<CalcHistIntegrals>("ihistogram");
//        Calculators().Register<ReadXML>("readxml");

	Calculators().Register<CalcEstatics>("ecoulomb");
        Calculators().Register<GenerateNrgs>("egaussian");
        Calculators().Register<EnergyCorr>("ecorrelation");
        Calculators().Register<CalcHistEnergeticDisorder>("ehistogram");
        Calculators().Register<Eshuffle>("eshuffle");
        Calculators().Register<Etinker>("etinker");
	
        Calculators().Register<Neighborlist>("neighborlist");
        Calculators().Register<Oboltzmann>("oboltzmann");
	Calculators().Register<Vaverage>("vaverage"); // average charge velocities

	Calculators().Register<CalcLambdaOut>("lambdaout");
        Calculators().Register<MarcusRates>("marcusrates");
        Calculators().Register<PolymerRates>("polymerrates");
        Calculators().Register<PairDump>("pairdump");
        Calculators().Register<CurrentDensity>("currden");
        Calculators().Register<JortnerRates>("jortnerrates");
        Calculators().Register<MarcusRatesLambdaOuter>("marcusrateslambdaouter");
        Calculators().Register<SQLiteWriter>("sqlitewriter");
        Calculators().Register<DumpAtomsBJ>("dumpatomsbj");
        Calculators().Register<DumpTrajectory>("tdump");
        Calculators().Register<Rates>("rates");
        Calculators().Register<WriteXML>("writexml");
	#ifdef WITH_VOTCA_KMCOLD        
        Calculators().Register<ContKmc>("kmc");
        #endif

}
