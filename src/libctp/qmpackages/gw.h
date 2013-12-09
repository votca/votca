/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __VOTCA_CTP_GW_H
#define	__VOTCA_CTP_GW_H


#include <votca/ctp/qmpackage.h>
#include <string> 

using namespace std;

namespace votca { namespace ctp {
/**
    \brief Wrapper for the Gaussian program
 
    The Gaussian class executes the Gaussian package 
    and extracts information from its log and io files
    
*/
class GW     : public QMPackage
{
public:   

   string getPackageName() { return "gw"; }

   void Initialize( Property *options );

   bool WriteInputFile( vector< Segment* > segments, Orbitals* _orbitals );

   bool Run();

   void CleanUp();
   
   bool CheckLogFile();

   bool ParseLogFile( Orbitals* _orbitals );

   bool ParseOrbitalsFile( Orbitals* _orbitals );

   bool ConvertToGW( Orbitals* _orbitals );
   
private:  

    string                              _scratch_dir;
    string                              _cleanup;
    string                              _ranges;
    string                              _gwbasis;
    
    double                              _rpamaxfactor;
    double                              _qpminfactor;
    double                              _qpmaxfactor;
    double                              _bseminfactor;
    double                              _bsemaxfactor;
    unsigned int                        _rpamax;
    unsigned int                        _qpmin;
    unsigned int                        _qpmax;
    unsigned int                        _bsemin;
    unsigned int                        _bsemax;

};


}}

#endif	/* __VOTCA_CTP_GW_H */

