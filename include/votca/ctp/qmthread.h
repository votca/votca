#ifndef _CTP_QMTHREAD_H
#define	_CTP_QMTHREAD_H

#include "votca/tools/thread.h"
#include "votca/tools/globals.h"
#include "votca/ctp/logger.h"
#include <iostream>
#include <string>
#include <ctime>

namespace votca { namespace ctp {
    

    // ++++++++++++++++++++++++++++++++++++++ //
    // Thread class with local string stream //
    // ++++++++++++++++++++++++++++++++++++++ //

     class QMThread : public Thread
    {

        /*
        friend ostream& operator<<( ostream& out, QMThread&  t ) {
           out << (t._ss).str();
           t._ss.str( "" );
           return out;
       }
       

        template <class Traits> 
        friend QMThread& operator<<( QMThread &t,  Traits& inp ) {
            
            if ( tools::globals::verbose ) { 
                if ( t._maverick ) { std::cout << inp ; }
                else { t._ss << inp ; }
            }
            return t;
        }
        
        template <class Traits> 
        friend QMThread& operator>>( Traits& inp, QMThread & t ) {
            
            if ( tools::globals::verbose ) { 
                if ( t._maverick ) { std::cout << inp ; }
                else { t._ss << inp ; }
            }           
            return t;
        }        
        */
       
    public:
        QMThread( bool maverick = false ) { _maverick = maverick; }; 
       ~QMThread() {};

       Logger* getLogger() { return &_logger; }
       
    protected:
        
        stringstream     _ss;
        bool             _maverick;
        Logger          _logger;

    };
    
    
}}

#endif /* _CTP_QMTHREAD_H */