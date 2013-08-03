/* 
 * File:   qmpackage.h
 * Author: andrienko
 *
 * Created on July 30, 2013, 8:25 PM
 */
#ifndef _CTP_QMPACKAGE_H
#define	_CTP_QMPACKAGE_H

#include <votca/ctp/logger.h>

namespace votca { namespace ctp {

using namespace std;    
// ========================================================================== //
// QMPackage base class for wrappers of TURBOMOLE, GAUSSIAN, etc          //
// ========================================================================== //
    
class QMPackage
{

public:

   QMPackage(){};
   virtual ~QMPackage(){};

   virtual string getPackageName() = 0;
   
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
    
private:
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
    
    Logger*                             _pLog;
    
};

}}

#endif	/* QMPACKAGE_H */

