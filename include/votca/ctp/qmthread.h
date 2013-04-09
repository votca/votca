#ifndef _CTP_QMTHREAD_H
#define	_CTP_QMTHREAD_H

#include "votca/tools/thread.h"
#include <iostream>
#include <string>

namespace votca { namespace ctp {

// ++++++++++++++++++++++++++++++++++++++ //
    // Thread class with local string stream //
    // ++++++++++++++++++++++++++++++++++++++ //

    class QMThread : public Thread
    {
       friend ostream& operator<<( ostream& output, const QMThread&  thread ) {
           output << (thread._ss).str();
           return output;
       }
       
       friend string operator>>( string input, QMThread&  thread ) {
           thread._ss << input ;
           return input;
       }
       
       friend int& operator>>( int& input, QMThread&  thread ) {
           thread._ss << input ;
           return input;
       }

       friend double& operator>>( double& input, QMThread&  thread ) {
           thread._ss << input ;
           return input;
       }

       friend float& operator>>( float& input, QMThread&  thread ) {
           thread._ss << input;
           return input;
       }
       
    public:
        
       ~QMThread() {};
       void AddLine( const string _line ) {
           _ss << _line;
       }

    protected:
        
        stringstream     _ss;
    };


    
    
}}

#endif 