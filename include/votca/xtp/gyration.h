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
#include <votca/ctp/logger.h>
#include <votca/xtp/qmmachine.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class Density2Gyration 
{
public:

    Density2Gyration (ctp::Logger* log) {_log=log; }
   ~Density2Gyration () { 
   
    std::vector< QMAtom* >::iterator it;
    for ( it = _Atomlist.begin(); it != _Atomlist.end(); ++it ) delete *it;};

    string Identify() { return "density2gyration"; }

    void   Initialize(Property *options);
    
    ub::vector<double> get_quaternion( ub::matrix<double> &eigenframe );
   
    void Convert2Eigenframe( ub::vector<double> V, ub::vector<double> &_diagonal, ub::matrix<double> &_eigenframe  );
    void ReportAnalysis( string label, ub::vector<double> _tensor_elements, ub::vector<double> _tensor_diagonal, ub::matrix<double> _tensor_frame );
    
    
    void AnalyzeDensity( Orbitals& _orbitals );
    void AnalyzeGeometry( vector< QMAtom* > _atoms );

private:
    
    int         _state_no;  
    int         _openmp_threads;
    string      _state;
    string      _method;
    string      _spin;
    string      _integrationmethod;
    string      _gridsize;



    vector< QMAtom* > _Atomlist;
    
    ctp::Logger*      _log;
    
    

};



}}



#endif /* _MUSCET_XTP_GYRATION_H */

