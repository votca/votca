/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef _XTP_QMTHREAD_H
#define	_XTP_QMTHREAD_H

#include "votca/tools/thread.h"
#include "votca/tools/globals.h"
#include "votca/xtp/logger.h"
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
        QMThread() { _maverick = false; }
        QMThread(bool maverick) { _maverick = maverick; }; 
       ~QMThread() {};

        int  getId() { return _id; }
        void setId(int id) { _id = id; }
        bool isMaverick() { return _maverick; }
       
        Logger* getLogger() { return &_logger; }
        virtual void Run(void) { ; }
       
    protected:
        
        int              _id;
        std::stringstream     _ss;
        bool             _maverick;
        Logger           _logger;

    };
    
    
}}

#endif /* _XTP_QMTHREAD_H */
