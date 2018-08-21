/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_ESP2MULTIPOLE_H
#define _VOTCA_XTP_ESP2MULTIPOLE_H

#include <stdio.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/mulliken.h>
#include <votca/xtp/lowdin.h>
#include <votca/xtp/nbo.h>
#include <votca/ctp/logger.h>
#include <boost/filesystem.hpp>
#include <votca/tools/property.h>
#include <votca/xtp/orbitals.h>

namespace votca { namespace xtp {

    
class Esp2multipole 
{
public:

    Esp2multipole (ctp::Logger* log) {
        _log=log; 
        _pairconstraint.resize(0);
        _regionconstraint.resize(0);
        }
   ~Esp2multipole () {};

    std::string Identify() { return "esp2multipole"; }

    void   Initialize(tools::Property &options);
    
   
    void Extractingcharges( Orbitals& orbitals );
    void WritetoFile(std::string output_file,  std::string identifier="esp2multipole");
    std::string GetIdentifier();

private:
    
    int         _state_no;  
    int         _openmp_threads;
    std::string      _state;
    std::string      _method;
    std::string      _spin;
    std::string      _integrationmethod;
    std::string      _gridsize;
    bool        _use_mulliken;
    bool        _use_lowdin;
    bool        _use_CHELPG;
    bool        _use_bulkESP;
    bool        _use_GDMA;
    bool        _use_CHELPG_SVD;
    bool        _use_NBO;
    bool        _use_ecp;
    bool        _do_svd;
    double      _conditionnumber;
    std::vector< QMAtom* > _Atomlist;
    
    ctp::Logger*      _log;
    std::vector< std::pair<int,int> > _pairconstraint; //  pairconstraint[i] is all the atomindices which have the same charge     
    std::vector< Espfit::region > _regionconstraint; 
    

};



}}


#endif