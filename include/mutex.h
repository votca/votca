/* 
 * File:   mutex.h
 * Author: koschke
 *
 * Created on October 28, 2010, 3:26 PM
 */

#ifndef MUTEX_H
#define	MUTEX_H
#include <pthread.h>

class Mutex {
public:
    Mutex();
    ~Mutex();

    //do we need recursive locks? i.e. locking a  mutex 3 times requires the same mutex to be unlocked
    //3 times

    void Lock();
    void Unlock();

private:
    pthread_mutex_t _mutexVar;
};


#endif	/* MUTEX_H */

