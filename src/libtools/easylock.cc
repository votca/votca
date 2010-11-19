/*
 * easylock.cc
 *
 *  Created on: Oct 26, 2010
 *      Author: koschke
 */

#include "easylock.h"

namespace votca { namespace tools {

EasyLock::EasyLock(Mutex *mutex) {
  // TODO Auto-generated constructor stub
  myMutex = mutex;
  myMutex->Lock();
}

EasyLock::~EasyLock() {
  // TODO Auto-generated destructor stub
  myMutex->Unlock();
}

}}