/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef _XTP_QM_PACKAGE_H
#define	_XTP_QM_PACKAGE_H

#include <votca/ctp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/tools/property.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/qmpair.h>
#include <votca/ctp/topology.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {
  
 
// ========================================================================== //
// QMPackage base class for wrappers of TURBOMOLE, GAUSSIAN, etc              //
// ========================================================================== //
    
class QMPackage
{

public:

   QMPackage(){};
   virtual ~QMPackage(){}; 

   virtual std::string getPackageName() = 0;


   virtual void Initialize( Property *options ) = 0;
   
   /// writes a coordinate file WITHOUT taking into account PBCs
   virtual bool WriteInputFile( std::vector< ctp::Segment* > segments, Orbitals* orbitals = NULL) = 0;

   /// writes a coordinate file of a pair WITH PBCs and the orbital guess [if needed]
   bool WriteInputFilePBC( ctp::QMPair* pair, Orbitals* orbitals = NULL);
   
   virtual bool Run() = 0;

   virtual bool ParseLogFile( Orbitals* _orbitals ) = 0;

   virtual bool ParseOrbitalsFile( Orbitals* _orbitals ) = 0;
   
   virtual void CleanUp() = 0;
   

   void setRunDir( std::string run_dir ) { _run_dir = run_dir; }
   
   void setInputFileName( std::string input_file_name ) { _input_file_name = input_file_name; }

   void setLogFileName( std::string log_file_name ) { _log_file_name = log_file_name; }

   void setOrbitalsFileName( string orb_file ) { _orb_file_name = orb_file; }
   
   void setLog( ctp::Logger* pLog ) { _pLog = pLog; }
      
   bool GuessRequested( ) { return _write_guess; }
   
   bool ECPRequested( ) { return _write_pseudopotentials; }
   
   bool VXCRequested() { return _output_Vxc; }

   void setCharge(const int charge) { _charge = charge; }
   
   void setSpin(const int spin) { _spin = spin; }
   
   void setThreads(const int threads) { _threads = threads; }
   
   void doGetCharges(bool do_get_charges) { _get_charges = do_get_charges; }
   
   std::string getBasisSetName(){return _basisset_name;}
   std::string getExecutable() {return _executable;};
   
protected:

    int                                 _charge;
    int                                 _spin; // 2S+1
    int                                 _threads;
    std::string                              _memory;
    std::string                              _options;
    
    std::string                              _executable;
    std::string                              _input_file_name;
    std::string                              _log_file_name;
    std::string                              _xyz_file_name;
    std::string                              _orb_file_name;
    
    std::string                              _run_dir;
        
    std::string                              _basisset_name;
    std::list< std::string >                      _cleanup_list;
    
    bool                                _get_orbitals;
    bool                                _get_overlap;
    bool                                _get_charges;
    bool                                _get_self_energy;
    
    bool                                _write_guess;
    bool                                _write_charges;
    bool                                _write_basis_set;
    bool                                _write_pseudopotentials;
    
    bool                                _output_Vxc;
    
    ctp::Logger*                             _pLog;
       
};

inline bool QMPackage::WriteInputFilePBC( ctp::QMPair* pair, Orbitals* orbitals) {
    
    //std::cout << "IDFT writes input with PBC" << std::endl;
    
    ctp::Segment* seg1 = pair->Seg1();
    ctp::Segment* seg2 = pair->Seg2();
    ctp::Segment* ghost = NULL;
    
    ctp::Topology* _top = seg1->getTopology();

    ctp::vec r1 = seg1->getPos();
    ctp::vec r2 = seg2->getPos();

    ctp::vec _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

    // Check whether pair formed across periodic boundary
    if ( abs(r2 - r1 - _R) > 1e-8 ) {
        ghost = new ctp::Segment(seg2);
        //ghost->TranslateBy(r1 - r2 + _R); // DO NOT USE THIS METHOD !
	std::vector<ctp::Atom*>::iterator ait;
	for (ait = ghost->Atoms().begin(); ait != ghost->Atoms().end(); ++ait) {
		(*ait)->setQMPos((*ait)->getQMPos()+r1-r2+_R);
	}
    }
 
    std::vector< ctp::Segment* > segments;
    segments.push_back(seg1);
    
    if ( ghost ) {
        segments.push_back(ghost);
    } else {
        segments.push_back(seg2);
    }
   
    WriteInputFile( segments, orbitals);
    
    delete ghost;
    return true;
}   

}}

#endif	/* _XTP_QM_PACKAGE_H */

