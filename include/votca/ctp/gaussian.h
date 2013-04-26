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

#ifndef __VOTCA_CTP_GAUSSIAN_H
#define	__VOTCA_CTP_GAUSSIAN_H

#include <votca/tools/property.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/apolarsite.h>

#include <string> 

using namespace std;

namespace votca { namespace ctp {
/**
    \brief Wrapper for the Gaussian program
 
    The Gaussian class executes the Gaussian package 
    and extracts information from its log and io files
    
*/
class Gaussian
{
public:   

    Gaussian( tools::Property *opt );
   ~Gaussian();

   /* Writes Gaussian input file with coordinates taken from all segments
    * and guess for the dimer orbitals (if given) constructed from the
    * orbitals of monomers 
    */
   bool WriteInputFile( vector< Segment* > segments, Orbitals* orbitals_guess = NULL);
   
   bool WriteShellScript();
   bool Run();
   void CleanUp( string ID );
   
   bool ParseLogFile( Orbitals* _orbitals );
   bool ParseOrbitalsFile( Orbitals* _orbitals );
   
   void setScratchDir( string scratch_dir ) { _scratch_dir = scratch_dir; }
   void setRunDir( string run_dir ) { _run_dir = run_dir; }
   void setInputFile( string com_file ) { _com_file_name = com_file; }
   void setShellFile( string shell_file ) { _shell_file_name = shell_file; }
   void setLogFile( string log_file ) { _log_file_name = log_file; }
   void setOrbitalsFile( string orb_file ) { _orb_file_name = orb_file; }
   
   string getScratchDir( ) { return _scratch_dir; }
   bool GuessRequested( ) { return _write_guess; }
   
//protected:
             
private:  

    static const double _conv_Hrt_eV = 27.21138386;
    
    int                                 _charge;
    int                                 _spin; // 2S+1
    string                              _options;
    
    string                              _executable;
    string                              _memory;
    int                                 _threads;
    
    string                              _shell_file_name;
    string                              _com_file_name;
    string                              _log_file_name;
    string                              _xyz_file_name;
    string                              _chk_file_name;
    string                              _orb_file_name;
    
    string                              _run_dir;
    string                              _scratch_dir;
        
    string                              _cleanup;
    
    bool                                _get_orbitals;
    bool                                _get_overlap;
    bool                                _get_charges;
    bool                                _get_self_energy;
    
    bool                                _write_guess;
    bool                                _write_charges;
         
    int NumberOfElectrons( string _line ); 
    int BasisSetSize( string _line ); 
    int EnergiesFromLog( string _line, ifstream inputfile ); 
    string FortranFormat( const double &number );
    
};


}}

#endif	/* __VOTCA_CTP_GAUSSAIN_H */

