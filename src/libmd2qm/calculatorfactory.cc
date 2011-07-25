#include <votca/ctp/calculatorfactory.h>

#include "votca_config.h"
#include "calculators/ecoulomb.h"
#include "calculators/ehistogram.h"
#include "calculators/egaussian.h"
#include "calculators/ecorrelation.h"
#include "calculators/eshuffle.h"
#include "calculators/eoutersphere.h"

#include "calculators/izindo.h"
#include "calculators/ihistogram.h"

#include "calculators/tdump.h"

#include "calculators/vaverage.h"

#include "calculators/oboltzmann.h"

#include "calculators/neighborlist.h"
#include "calculators/etinker.h"

#include "calculators/rates.h"

#include "calculators/writexml.h"
#include "calculators/polymerrates.h"
#include "calculators/pairdump.h"
#include "calculators/CurrentDensity.h"
#include "calculators/sqlitewriter.h"

//#include "calculators/readxml.h"
//#include "calculators/marcusrates.h"
//#include "calculators/jortnerrates.h"
//#include "calculators/marcusrateslambdaouter.h"

#ifdef WITH_VOTCA_KMCOLD        
    #include "calculators/contkmc.h"
#endif

void CalculatorFactory::RegisterAll(void)
{
	Calculators().Register<CalcIntegrals>("izindo"); // ZINDO-based transfer integrals
        Calculators().Register<CalcHistIntegrals>("ihistogram"); // histogram of transfer integrals

	Calculators().Register<CalcEstatics>("ecoulomb"); // Coulomb part of site energies
	Calculators().Register<CalcLambdaOut>("eoutersphere"); // Outersphere reorganization energy
        Calculators().Register<Eshuffle>("eshuffle"); // removes spatial energy correlations
        Calculators().Register<GenerateNrgs>("egaussian"); // gaussian distribution of site energies
        Calculators().Register<EnergyCorr>("ecorrelation"); // site energy correlation function
        Calculators().Register<CalcHistEnergeticDisorder>("ehistogram"); // site energy histogram
        Calculators().Register<Etinker>("etinker"); // input for the TINKER package (site energies)
	
        Calculators().Register<Neighborlist>("neighborlist"); // fragment-based neighbor list
        Calculators().Register<Oboltzmann>("oboltzmann"); // Boltzmann distribution of site energies
	Calculators().Register<Vaverage>("vaverage"); // average charge velocities (requires site occupations)

        Calculators().Register<PolymerRates>("polymerrates");
        Calculators().Register<PairDump>("pairdump");
        Calculators().Register<CurrentDensity>("currden");
        Calculators().Register<SQLiteWriter>("sqlitewriter");

	Calculators().Register<DumpTrajectory>("tdump"); // coarse-grained and based on rigid segments trajectories
        Calculators().Register<Rates>("rates"); // Marcus, Jortner rates
        Calculators().Register<WriteXML>("writexml");  // obsolete

	#ifdef WITH_VOTCA_KMCOLD        
        Calculators().Register<ContKmc>("kmc");
        #endif

}
