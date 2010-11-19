/* 
 * File:   thread.h
 * Author: koschke
 *
 * Created on October 28, 2010, 3:11 PM
 */

#ifndef THREAD_H
#define	THREAD_H

#include <pthread.h>

namespace votca { namespace tools {

class Thread {
public:
    Thread();
    virtual ~Thread();

    void Start();
    virtual void Run(void) = 0;

    void WaitDone();
    bool IsFinished() const;

    //wait for others
    //syncronization: using mutex (except if the order of the threads is important - use signals/conditions)

    //signals
    //void Finished();

    //function that automatically determines the number of threads (number of cpus? # cpus * 2, if HT or SMT is available?)
    //static int getNumberOfThreadsToUse();

private:
    pthread_t _thread;
    bool _finished;

};

}}

#endif	/* THREAD_H */

