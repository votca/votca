/*
 * easylock.h
 *
 *  Created on: Oct 26, 2010
 *      Author: koschke
 */

#ifndef EASYLOCK_H_
#define EASYLOCK_H_

#include "mutex.h"

namespace votca { namespace tools {

class EasyLock {
public:
  EasyLock(Mutex *mutex);
  //destructor unlocks mutex
  ~EasyLock();

  //return mutex used for lock
  //Mutex * mutex() const;

  //void Unlock();
  //void Relock();

private:
  Mutex *myMutex;
};

}}

#endif /* EASYLOCK_H_ */
