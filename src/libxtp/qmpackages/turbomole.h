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

#ifndef __VOTCA_XTP_TURBOMOLE_H
#define	__VOTCA_XTP_TURBOMOLE_H

#include <votca/tools/property.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/apolarsite.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmpackage.h>

#include <string> 

using namespace std;

namespace votca { namespace xtp {
/**
    \brief Wrapper for the Gaussian program
 
    The Gaussian class executes the Gaussian package 
    and extracts information from its log and io files
    
*/
class Turbomole : public QMPackage
{
public:   

   string getPackageName() { return "turbomole"; }

   void Initialize( Property *options );

   /* Writes Turbomole input file with coordinates taken from all segments
    * and guess for the dimer orbitals (if given) constructed from the
    * orbitals of monomers 
    */
   bool WriteInputFile( vector< Segment* > segments, Orbitals* orbitals_guess = NULL);

   bool Run();

   void CleanUp();
   
   bool CheckLogFile();

   bool ParseLogFile( Orbitals* _orbitals );

   bool ParseOrbitalsFile( Orbitals* _orbitals );

   bool ConvertToGW( Orbitals* _orbitals );
   
private:  

    //static const double _conv_Hrt_eV = 27.21138386;
    
    //int                                 _charge;
    //int                                 _spin; // 2S+1
    string                              _options;
    
    string                              _memory;
    //int                                 _threads;
    
    string                              _scratch_dir;
        
    string                              _cleanup;

    string FortranFormat( const double &number );    
};


}}

#endif	/* __VOTCA_XTP_GAUSSAIN_H */

