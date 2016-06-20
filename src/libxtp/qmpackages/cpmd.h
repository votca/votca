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

#ifndef __VOTCA_XTP_CPMD_H
#define	__VOTCA_XTP_CPMD_H


#include <votca/xtp/apolarsite.h>
#include <votca/xtp/qmpackage.h>

#include <string> 



namespace votca { namespace xtp {
/**
    \brief Wrapper for the CPMD program
 
    The Cpmd class executes the CPMD package 
    and extracts information from its log and io files
    
*/
class Cpmd : public QMPackage
{
public:   

   std::string getPackageName() { return "cpmd"; }

   void Initialize( Property *options );

   /* Writes CPMD input file */
   bool WriteInputFile( std::vector< Segment* > segments, Orbitals* orbitals_guess = NULL);

   bool Run();

   void CleanUp();
   
   bool CheckLogFile();

   bool ParseLogFile( Orbitals* _orbitals );
   
private:  

    std::string                              _shell_file_name;
    std::string                              _chk_file_name;
    std::string                              _scratch_dir;
    std::string                              _input_vxc_file_name;    
    std::string                              _cleanup;
    std::string                              _vdWfooter;

    int NumberOfElectrons( std::string _line ); 
    int BasisSetSize( std::string _line ); 
    int EnergiesFromLog( std::string _line, std::ifstream inputfile ); 
    std::string FortranFormat( const double &number );
    int NumbfQC( std::string _shell_type);
    int NumbfGW( std::string _shell_type);
    int NumbfQC_cart( std::string _shell_type);

    
    
};


}}

#endif	/* __VOTCA_XTP_CPMD_H */

