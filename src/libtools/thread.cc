/* 
 * File:   thread.cc
 * Author: koschke
 * 
 * Created on October 28, 2010, 3:11 PM
 */

#include "thread.h"
#include <stdexcept>

static void *runwrapper(void *arg) {
    Thread *thread = (Thread*) (arg);
    thread->Run();
    pthread_exit(NULL);
}

Thread::Thread() {

}

Thread::~Thread() {
    // TODO Auto-generated destructor stub

}

void Thread::Start() {
    _finished = false;

    int rc = pthread_create(&_thread, NULL, runwrapper, (void *) this);
    if (rc) {
        throw std::runtime_error("ERROR; return code from pthread_create() is "
                + rc);
    }
    //WaitDone();
    //_finished=true;
}

void Thread::WaitDone() {
    //is there a better implementation using conditions?
    void * status;
    int rc = pthread_join(_thread, &status);
    if (rc) {
        throw std::runtime_error("ERROR; return code from pthread_join() is "
                + rc);
    }
    _finished = true;
}

bool Thread::IsFinished() const {
    return _finished;
    //	if (_finished)
    //		return true;
    //	else
    //		return false;
}

