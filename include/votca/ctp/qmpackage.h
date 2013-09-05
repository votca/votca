/* 
 * File:   qmpackage.h
 * Author: andrienko
 *
 * Created on July 30, 2013, 8:25 PM
 */
#ifndef _CTP_QM_PACKAGE_H
#define	_CTP_QM_PACKAGE_H

#include <votca/ctp/logger.h>
#include <votca/ctp/orbitals.h>
#include <votca/tools/property.h>
#include <votca/ctp/segment.h>

namespace votca { namespace ctp {

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
   
   virtual bool WriteInputFile( vector< Segment* > segments, Orbitals* orbitals = NULL) = 0;
   
   virtual bool Run() = 0;

   virtual bool ParseLogFile( Orbitals* _orbitals ) = 0;

   virtual bool ParseOrbitalsFile( Orbitals* _orbitals ) = 0;
   
   virtual void CleanUp() = 0;

   void setRunDir( string run_dir ) { _run_dir = run_dir; }
   
   void setInputFileName( string input_file_name ) { _input_file_name = input_file_name; }

   void setLogFileName( string log_file_name ) { _log_file_name = log_file_name; }

   void setOrbitalsFileName( string orb_file ) { _orb_file_name = orb_file; }
   
   void setLog( Logger* pLog ) { _pLog = pLog; }
      
   bool GuessRequested( ) { return _write_guess; }

   void setCharge(const int charge) { _charge = charge; }
   
   void setSpin(const int spin) { _spin = spin; }
   
   void setThreads(const int threads) { _threads = threads; }
   
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
        
    list< string >                      _cleanup_list;
    
    bool                                _get_orbitals;
    bool                                _get_overlap;
    bool                                _get_charges;
    bool                                _get_self_energy;
    
    bool                                _write_guess;
    bool                                _write_charges;
    bool                                _is_optimization;
    
    Logger*                             _pLog;
       
};

}}

#endif	/* _CTP_QM_PACKAGE_H */

