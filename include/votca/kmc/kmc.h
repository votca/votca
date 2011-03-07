#ifndef __VOTCA_KMC_KMC_H
#define __VOTCA_KMC_KMC_H

#include <iostream>
#include <queue>

namespace votca { namespace kmc {

using namespace std;

/**
  \brief implements first reaction method kmc
*/
template<typename event_type>
class Kmc
{
public:

    class event {
        event(double time, int id)
            : _time(time), _id(id) {};

        void execute(Kmc &kmc) {
            cout << "event: " << _id 
                 << " time: " << _time 
                 << endl;
        }

    protected:
        // time when this event occurs
        double _time;
        // identifier for the event
        int _id;
    };

    double getSystemTime() { return _system_time; }
    void addEvent(event_type &event);

protected:
    double _system_time;

    // all events are stored in a priority queue
    priority_queue<event_type> _events;
};


template<typename event_type>
void Kmc<event_type>::addEvent(event_type &event)
{
    _events.push(event);
}

}}
#endif

