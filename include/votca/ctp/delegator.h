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

#ifndef __VOTCA_CTP_DELEGATOR_H
#define	__VOTCA_CTP_DELEGATOR_H

#include <iostream>
#include <boost/interprocess/sync/file_lock.hpp> 
using namespace boost::interprocess ; 

namespace votca { namespace ctp {


enum TJobStatus {jobAVAILABLE, jobLOCKED, jobDONE, jobFAILED};

 
/**
*   \brief Delegator distributes jobs
*/
class Delegator { 
    
public:
        Delegator() {}
       ~Delegator() {}
        
        void updateJobList() {
            std::ofstream os ;
            os.open(m_FileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::app) ; 
            file_lock flock("lock");
            if (flock.try_lock()) { // flock.lock() - with waiting // timed_lock
                cout << "Acquired lock" << endl; }
            else {
                cout << "Could not acquire lock" << endl; 
                return 1;
            }
            sleep(30);
            os.close();
        }   
        void getJobs() {}
        void giveJob(){}
        
private:
        std:map _map<int,string>;
        std::string m_FileName("./locktest.txt") ; 
        
};

}}

#endif /* __VOTCA_CTP_DELEGATOR_H */