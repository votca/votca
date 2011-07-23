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

#include "calculators/rates.h"

//#include "calculators/readxml.h"
//#include "calculators/marcusrates.h"
//#include "calculators/jortnerrates.h"
//#include "calculators/marcusrateslambdaouter.h"

#include "calculators/writexml.h"
#include "calculators/polymerrates.h"
#include "calculators/pairdump.h"
#include "calculators/CurrentDensity.h"
#include "calculators/sqlitewriter.h"
#include "calculators/lambdaout.h"
#include "calculators/dump_atoms_bj.h"

#ifdef WITH_VOTCA_KMCOLD        
    #include "calculators/contkmc.h"
#endif

void CalculatorFactory::RegisterAll(void)
{
	Calculators().Register<CalcIntegrals>("izindo"); // ZINDO-based transfer integrals
        Calculators().Register<CalcHistIntegrals>("ihistogram"); // histogram of transfer integrals

	Calculators().Register<CalcEstatics>("ecoulomb"); // Coulomb part of site energies
        Calculators().Register<GenerateNrgs>("egaussian"); // gaussian distribution of site energies
        Calculators().Register<EnergyCorr>("ecorrelation"); // site energy correlation function
        Calculators().Register<CalcHistEnergeticDisorder>("ehistogram"); // site energy histogram
        Calculators().Register<Eshuffle>("eshuffle"); // removes spatial energy correlations
        Calculators().Register<Etinker>("etinker"); // input for the TINKER package (site energies)
	
        Calculators().Register<Neighborlist>("neighborlist"); // fragment-based neighbor list
        Calculators().Register<Oboltzmann>("oboltzmann"); // Boltzmann distribution of site energies
	Calculators().Register<Vaverage>("vaverage"); // average charge velocities (requires site occupations)

	Calculators().Register<CalcLambdaOut>("lambdaout");
//        Calculators().Register<MarcusRates>("marcusrates");
        Calculators().Register<PolymerRates>("polymerrates");
        Calculators().Register<PairDump>("pairdump");
        Calculators().Register<CurrentDensity>("currden");
//        Calculators().Register<JortnerRates>("jortnerrates");
//        Calculators().Register<MarcusRatesLambdaOuter>("marcusrateslambdaouter");
        Calculators().Register<SQLiteWriter>("sqlitewriter");
        Calculators().Register<DumpAtomsBJ>("dumpatomsbj");
        Calculators().Register<DumpTrajectory>("tdump");
        Calculators().Register<Rates>("rates");
        Calculators().Register<WriteXML>("writexml");
//        Calculators().Register<ReadXML>("readxml");
	#ifdef WITH_VOTCA_KMCOLD        
        Calculators().Register<ContKmc>("kmc");
        #endif

}
