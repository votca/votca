/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#include <votca/ctp/toolfactory.h>
#include <votca/xtp/toolfactory.h>
#include <votca/xtp/qmtool.h>



#include "tools/molpol.h"
/*#include "tools/pdb2map.h"
#include "tools/coupling.h"
#include "tools/log2mps.h"
#include "tools/ptopreader.h"
#include "tools/pdb2top.h"
#include "tools/exciton.h"
#include "tools/qmanalyze.h"
#include "tools/qmsandbox.h"
#include "tools/spectrum.h"
#include "tools/excitoncoupling.h"
#include "tools/orb2isogwa.h"
#include "tools/dft.h"
#include "tools/gencube.h"
#include "tools/partialcharges.h"
*/

namespace XTP = votca::xtp;
namespace CTP = votca::ctp;

// votca::xtp::XQMToolFactory::

//namespace votca { namespace xtp {

void votca::xtp::QMToolFactory::RegisterAll(void)
{
        XTP::QMTools().Register<XTP::MolPolTool>         ("molpol");
        /*QMTools().Register<XTP::PDB2Map>            ("pdb2map");
        QMTools().Register<XTP::Coupling>           ("coupling");
        QMTools().Register<XTP::Log2Mps>            ("log2mps");
        QMTools().Register<XTP::PtopReader>         ("ptopreader");
        QMTools().Register<XTP::Exciton>            ("exciton");
        QMTools().Register<XTP::QMAnalyze>          ("qmanalyze");
        QMTools().Register<XTP::QMSandbox>          ("qmsandbox");
        QMTools().Register<XTP::Spectrum>           ("spectrum");
        QMTools().Register<XTP::ExcitonCoupling>    ("excitoncoupling");
        QMTools().Register<XTP::Orb2IsoGWA>         ("orb2isogwa"); 
        QMTools().Register<XTP::PDB2Top>            ("pdb2top");
        QMTools().Register<XTP::DFT>                ("dft");
        QMTools().Register<XTP::GenCube>            ("gencube");
        QMTools().Register<XTP::Partialcharges>     ("partialcharges");
        */
}

//}}
