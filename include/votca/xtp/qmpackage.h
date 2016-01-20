/* 
 * File:   qmpackage.h
 * Author: andrienko
 *
 * Created on July 30, 2013, 8:25 PM
 */
#ifndef _XTP_QM_PACKAGE_H
#define	_XTP_QM_PACKAGE_H

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/tools/property.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/topology.h>

namespace votca { namespace xtp {

using namespace std;    
// ========================================================================== //
// QMPackage base class for wrappers of TURBOMOLE, GAUSSIAN, etc              //
// ========================================================================== //
    
class QMPackage
{

public:

   QMPackage(){};
   virtual ~QMPackage(){}; 

   virtual string getPackageName() = 0;


   virtual void Initialize( Property *options ) = 0;
   
   /// writes a coordinate file WITHOUT taking into account PBCs
   virtual bool WriteInputFile( vector< Segment* > segments, Orbitals* orbitals = NULL) = 0;

   /// writes a coordinate file of a pair WITH PBCs and the orbital guess [if needed]
   bool WriteInputFilePBC( QMPair* pair, Orbitals* orbitals = NULL);
   
   virtual bool Run() = 0;

   virtual bool ParseLogFile( Orbitals* _orbitals ) = 0;

   virtual bool ParseOrbitalsFile( Orbitals* _orbitals ) = 0;
   
   virtual void CleanUp() = 0;
   
   virtual bool ConvertToGW( Orbitals* _orbitals ) = 0;

   void setRunDir( string run_dir ) { _run_dir = run_dir; }
   
   void setInputFileName( string input_file_name ) { _input_file_name = input_file_name; }

   void setLogFileName( string log_file_name ) { _log_file_name = log_file_name; }

   void setOrbitalsFileName( string orb_file ) { _orb_file_name = orb_file; }
   
   void setLog( Logger* pLog ) { _pLog = pLog; }
      
   bool GuessRequested( ) { return _write_guess; }
   
   bool ECPRequested( ) { return _write_pseudopotentials; }
   
   bool VXCRequested() { return _output_Vxc; }

   void setCharge(const int charge) { _charge = charge; }
   
   void setSpin(const int spin) { _spin = spin; }
   
   void setThreads(const int threads) { _threads = threads; }
   
   void doGetCharges(bool do_get_charges) { _get_charges = do_get_charges; }
   
   string getBasisSetName(){return _basisset_name;}
   string getExecutable() {return _executable;};
   
protected:

    int                                 _charge;
    int                                 _spin; // 2S+1
    int                                 _threads;
    string                              _memory;
    string                              _options;
    
    string                              _executable;
    string                              _input_file_name;
    string                              _log_file_name;
    string                              _xyz_file_name;
    string                              _orb_file_name;
    
    string                              _run_dir;
        
    string                              _basisset_name;
    list< string >                      _cleanup_list;
    
    bool                                _get_orbitals;
    bool                                _get_overlap;
    bool                                _get_charges;
    bool                                _get_self_energy;
    
    bool                                _write_guess;
    bool                                _write_charges;
    bool                                _write_basis_set;
    bool                                _write_pseudopotentials;
    
    bool                                _output_Vxc;
    
    Logger*                             _pLog;
       
};

inline bool QMPackage::WriteInputFilePBC( QMPair* pair, Orbitals* orbitals) {
    
    //std::cout << "IDFT writes input with PBC" << std::endl;
    
    Segment* seg1 = pair->Seg1();
    Segment* seg2 = pair->Seg2();
    Segment* ghost = NULL;
    
    Topology* _top = seg1->getTopology();

    vec r1 = seg1->getPos();
    vec r2 = seg2->getPos();

    vec _R = _top->PbShortestConnect(r1, r2); // => _R points from 1 to 2

    // Check whether pair formed across periodic boundary
    if ( abs(r2 - r1 - _R) > 1e-8 ) {
        ghost = new Segment(seg2);
        //ghost->TranslateBy(r1 - r2 + _R); // DO NOT USE THIS METHOD !
	vector<Atom*>::iterator ait;
	for (ait = ghost->Atoms().begin(); ait != ghost->Atoms().end(); ++ait) {
		(*ait)->setQMPos((*ait)->getQMPos()+r1-r2+_R);
	}
    }
 
    vector< Segment* > segments;
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

