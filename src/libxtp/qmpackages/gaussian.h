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

#ifndef __VOTCA_XTP_GAUSSIAN_H
#define	__VOTCA_XTP_GAUSSIAN_H


#include <votca/xtp/apolarsite.h>
#include <votca/xtp/qmpackage.h>

#include <string> 

using namespace std;

namespace votca { namespace xtp {
/**
    \brief Wrapper for the Gaussian program
 
    The Gaussian class executes the Gaussian package 
    and extracts information from its log and io files
    
*/
class Gaussian : public QMPackage
{
public:   

   string getPackageName() { return "gaussian"; }

   void Initialize( Property *options );

   /* Writes Gaussian input file with coordinates of segments
    * and a guess for the dimer (if requested) constructed from the
    * monomer orbitals
    */
   bool WriteInputFile( vector< Segment* > segments, Orbitals* orbitals_guess = NULL);

   bool WriteShellScript();

   bool Run();

   void CleanUp();
   
   bool CheckLogFile();

   bool ParseLogFile( Orbitals* _orbitals );

   bool ParseOrbitalsFile( Orbitals* _orbitals );
   
   bool ConvertToGW( Orbitals* _orbitals );
      
   string getScratchDir( ) { return _scratch_dir; }
   
private:  

    string                              _shell_file_name;
    string                              _chk_file_name;
    string                              _scratch_dir;
    string                              _input_vxc_file_name;    
    string                              _cleanup;
    string                              _vdWfooter;

    int NumberOfElectrons( string _line ); 
    int BasisSetSize( string _line ); 
    int EnergiesFromLog( string _line, ifstream inputfile ); 
    string FortranFormat( const double &number );
    int NumbfQC( string _shell_type);
    int NumbfGW( string _shell_type);
    int NumbfQC_cart( string _shell_type);

    
    
};


}}

#endif	/* __VOTCA_XTP_GAUSSAIN_H */

