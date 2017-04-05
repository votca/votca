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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>


#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/xtp/elements.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/qminterface.h>

using boost::format;

namespace votca {
    namespace xtp {

    ctp::APolarSite *QMMInterface::Convert(ctp::QMAtom *atm, int id) {
        double A_to_nm = 0.1;
        vec pos = A_to_nm*vec(atm->x, atm->y, atm->z);
        double q = atm->charge;
        std::string elem = atm->type;
        double pol = 0.0;
        try {
            pol = _polar_table.at(elem);
        }
        catch(out_of_range) {
            std::cout << std::endl << "QMMInterface - no default polarizability given "
                << "for element type '" << elem << "'. Defaulting to 1A**3" << std::flush;
            pol = 1e-3;
        }

        ctp::APolarSite *new_aps = new ctp::APolarSite(id, elem);
        new_aps->setRank(0);
        new_aps->setPos(pos);
        new_aps->setQ00(q,0); // <- charge state 0 <> 'neutral'
        new_aps->setIsoP(pol);
        
        return new_aps;
    }
        
    ctp::PolarSeg *QMMInterface::Convert(std::vector<ctp::QMAtom*> &atms) {        
        ctp::PolarSeg *new_pseg = new ctp::PolarSeg();
        std::vector<ctp::QMAtom*>::iterator it;
        for (it = atms.begin(); it < atms.end(); ++it) {
            ctp::APolarSite *new_site = this->Convert(*it);
            new_pseg->push_back(new_site);
        }
        return new_pseg;
    }


    }
}
