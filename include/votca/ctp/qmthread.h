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
       friend ostream& operator<<( ostream& out, QMThread&  t ) {
           //out << "... ... ... ... ... ... ... ... ... ... ..." << endl;
           out << (t._ss).str();
           //out << "... ... ... ... ... ... ... ... ... ... ..." << endl;
           t._ss.str( "" );
           return out;
       }
       

        template <class Traits> 
        friend Traits& operator<<( QMThread & t,  Traits& inp ) {
            t._ss << inp ;
            return inp;
        }        
        
        template <class Traits> 
        friend Traits& operator>>( Traits& inp, QMThread & t ) {
            t._ss << inp ;
            return inp;
        }        
        
       
    public:
        
       ~QMThread() {};
       
    protected:
        
        stringstream     _ss;
    };


    
    
}}

#endif /* _CTP_QMTHREAD_H */