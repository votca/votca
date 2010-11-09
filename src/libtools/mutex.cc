/* 
 * File:   mutex.cc
 * Author: koschke
 * 
 * Created on October 28, 2010, 3:26 PM
 */

#include "mutex.h"

Mutex::Mutex() {
  // TODO Auto-generated constructor stub
  pthread_mutex_init(&_mutexVar, NULL);

}

Mutex::~Mutex() {
  // TODO Auto-generated destructor stub
  pthread_mutex_destroy(&_mutexVar);
}

void Mutex::Lock()
{
  pthread_mutex_lock(&_mutexVar);
}


void Mutex::Unlock()
{
  pthread_mutex_unlock(&_mutexVar);
}

