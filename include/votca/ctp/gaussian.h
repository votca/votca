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

#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/ctp/segment.h>
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

   bool WriteInputFile( Segment *seg );
   bool WriteShellScript();
   bool Run();
   
   void setScratchDir( string scratch ) { _scratch = scratch; }
   void setRunDir( string run_dir ) { _run_dir = run_dir; }
   void setInputFile( string com_file ) { _com_file_name = com_file; }
   void setShellFile( string shell_file ) { _shell_file_name = shell_file; }
   string getScratchDir( ) { return _scratch; }
   
//protected:
             
private:  
    
    string                              _functional;
    string                              _basis_set;
    int                                 _charge;
    int                                 _spin; // 2S+1
    string                              _options;
    
    string                              _executable;
    string                              _memory;
    int                                 _threads;
    string                              _checkpointfile;
    string                              _scratch;
    
    string                              _root_dir;
    string                              _source_file;
    
    string                              _shell_file_name;
    string                              _com_file_name;
    string                              _run_dir;
    
};


}}

#endif	/* __VOTCA_CTP_GAUSSAIN_H */

