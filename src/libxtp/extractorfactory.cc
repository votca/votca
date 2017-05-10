/*
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/extractorfactory.h>

#include "extractors/energyextractor.h"
#include "extractors/integralsextractor.h"
#include "extractors/ratesextractor.h"
#include "extractors/trajextractor.h"
#include "extractors/segmentsextractor.h"
#include "extractors/pairsextractor.h"
#include "extractors/occupationsextractor.h"



namespace votca { namespace xtp {

void ExtractorFactory::RegisterAll(void)
{	
        Extractors().Register<EnergyExtractor>             ("xenergy2xml");
        Extractors().Register<IntegralsExtractor>          ("xintegrals2xml");
        Extractors().Register<RatesExtractor>              ("xrates2xml");
        Extractors().Register<OccupationsExtractor>        ("xoccupations2xml");
        Extractors().Register<TrajExtractor>               ("xtrajectory2pdb");
        Extractors().Register<SegmentsExtractor>           ("xsegments2xml");
        Extractors().Register<PairsExtractor>              ("xpairs2xml");
}

}}
