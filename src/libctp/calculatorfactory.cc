/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/ctp/calculatorfactory.h>

#include "votca_config.h"

#include "calculators/izindo.h"
#include "calculators/ihistogram.h"

#include "calculators/ecoulomb.h"
#include "calculators/ehistogram.h"
#include "calculators/egaussian.h"
#include "calculators/ecorrelation.h"
#include "calculators/eoutersphere.h"
#include "calculators/etinker.h"

#include "calculators/tdump.h"
#include "calculators/pairdump.h"

#include "calculators/vaverage.h" 

#include "calculators/oboltzmann.h"

#include "calculators/neighborlist.h"

#include "calculators/rates.h"

//#include "calculators/writexml.h"

#include "calculators/sqlitewriter.h"

#ifdef WITH_VOTCA_KMCOLD        
    #include "calculators/contkmc.h"
#endif

namespace votca { namespace ctp {

void CalculatorFactory::RegisterAll(void)
{
	Calculators().Register<Izindo>("izindo"); // ZINDO-based transfer integrals
        Calculators().Register<CalcHistIntegrals>("ihistogram"); // histogram of transfer integrals

	Calculators().Register<Ecoulomb>("ecoulomb"); // Coulomb part of site energies
	Calculators().Register<Eoutersphere>("eoutersphere"); // Outersphere reorganization energy
	Calculators().Register<Egaussian>("egaussian"); // gaussian (also correlated) distribution of site energies
        Calculators().Register<Ecorrelation>("ecorrelation"); // site energy correlation function
        Calculators().Register<Ehistogram>("ehistogram"); // site energy histogram
        Calculators().Register<Etinker>("etinker"); // input for the TINKER package (site energies)
	
        Calculators().Register<Neighborlist>("neighborlist"); // fragment-based neighbor list
        Calculators().Register<Oboltzmann>("oboltzmann"); // Boltzmann distribution of site energies

        Calculators().Register<PairDump>("pairdump"); // Dumps pair coordinates for DFT coupling elements
	Calculators().Register<Vaverage>("vaverage"); // average charge velocities (via site occupations)
        Calculators().Register<SQLiteWriter>("sqlitewriter");

	Calculators().Register<Tdump>("tdump"); // coarse-grained and based on rigid segments trajectories
        Calculators().Register<Rates>("rates"); // Marcus, Jortner rates
        
        //Calculators().Register<WriteXML>("writexml");  // obsolete
 
	#ifdef WITH_VOTCA_KMCOLD        
        Calculators().Register<ContKmc>("kmc");
        #endif

}

}}