/* 
 *            Copyright 2016 The MUSCET Development Team
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


#ifndef _VOTCA_XTP_GYRATION_H
#define _VOTCA_XTP_GYRATION_H


#include <stdio.h>
#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>
#include <boost/filesystem.hpp>
#include <votca/xtp/numerical_integrations.h>
namespace votca { namespace xtp {
    
    
class Density2Gyration 
{
public:

    Density2Gyration (ctp::Logger* log) {_log=log; }

    std::string Identify() { return "density2gyration"; }

    void  Initialize( tools::Property &options);

    void AnalyzeDensity( Orbitals& orbitals );
    
private:
    void ReportAnalysis( std::string label,Gyrationtensor gyro, tools::matrix::eigensystem_t system );
    void AnalyzeGeometry( std::vector< QMAtom* > atoms );
      
    int         _openmp_threads;
    QMState      _state;
    bool _dostateonly;
    std::string      _integrationmethod;
    std::string      _gridsize;
    ctp::Logger*      _log;
    
    

};



}}



#endif /* _MUSCET_XTP_GYRATION_H */

