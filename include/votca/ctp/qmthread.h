#ifndef _CTP_QMTHREAD_H
#define	_CTP_QMTHREAD_H

#include "votca/tools/thread.h"
#include "votca/tools/globals.h"
#include <iostream>
#include <string>

namespace votca { namespace ctp {

    // ++++++++++++++++++++++++++++++++++++++ //
    // Thread class with local string stream //
    // ++++++++++++++++++++++++++++++++++++++ //

    class QMThread : public Thread
    {
       friend ostream& operator<<( ostream& out, QMThread&  t ) {
           out << (t._ss).str();
           t._ss.str( "" );
           return out;
       }
       

        template <class Traits> 
        friend Traits& operator<<( QMThread & t,  Traits& inp ) {
            
            if ( tools::globals::verbose ) { 
                if ( t._maverick ) { std::cout << inp ; }
                else { t._ss << inp ; }
            }
            return inp;
        }        
        
        template <class Traits> 
        friend Traits& operator>>( Traits& inp, QMThread & t ) {
            
            if ( tools::globals::verbose ) { 
                if ( t._maverick ) { std::cout << inp ; }
                else { t._ss << inp ; }
            }           
            return inp;
        }        
        
       
    public:
        QMThread( bool maverick = false ) { _maverick = maverick; };        
       ~QMThread() {};
       
    protected:
        
        stringstream     _ss;
        bool             _maverick;

    };


    
    
}}

#endif /* _CTP_QMTHREAD_H */