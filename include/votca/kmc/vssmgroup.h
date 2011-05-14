#ifndef __VOTCA_KMC_VSSMGROUP_H_
#define __VOTCA_KMC_VSSMGROUP_H_

#include <votca/tools/random.h>
#include <vector>

namespace votca { namespace kmc {

using namespace votca::tools;

/***
 * \brief groups events based on variable step size method
 *
 * This class groups the events based on the variable step size method. The waiting time is calculated
 * from the total rate. The probability to choose an event is P_i = \omega_i / \omega_tot
 *
 * The template argument is a class which has to implement the functions Rate, which should return the rate of the event
 * and onExecute, which is called if the event occurs. VSSMGroup implements these methods it's own and therefore
 * can be used as an event itself. This allows hierarchical grouping of events.
 *
 * This class is the most generic case and can handle changing rates on the fly. If the rates do not change
 * use VSSMStatic which has a better performance.
 *
 */
template<typename event_t>
class VSSMGroup {
public:
	VSSMGroup() {}

	/***
	 * \brief add event to list of possible events
	 */void AddEvent(event_t *event) {
		_events.push_back(event);
		onEventAdded(event);
		UpdateWaitingTime();
	}

	/**
	 * \brief get the total rate of this group
	 */
	double Rate();
	/**
	 * \brief this group is called, select and execute an event in group
	 */
	void onExecute();
	/**
	 * \brief get the waiting time for this group of events
	 */
	double WaitingTime();

	/**
	 * \brief calculate a new waiting time
	 */
	void UpdateWaitingTime() {
		_waiting_time = -log( 1.0 - Random::rand_uniform() ) / Rate();
	}

	event_t *getEvent(int i) { return _events[i]; }

protected:

	/**
	 * \brief this is called after AddEvent is finished
	 */
	void onEventAdded(event_t *event) {};

	/**
	 * \brief select an event to execute
	 */
	event_t *SelectEvent();

	std::vector<event_t*> _events;
	double _waiting_time;
};

template<typename event_t>
inline double VSSMGroup<event_t>::Rate()
{
	double r=0;
	for(typename std::vector<event_t*>::iterator e=_events.begin();
			e!=_events.end(); ++e) {
			r+=(*e)->Rate();
	}
	return r;
}

template<typename event_t>
inline void VSSMGroup<event_t>::onExecute()
{
	SelectEvent()->onExecute();
	UpdateWaitingTime();
}

template<typename event_t>
inline event_t *VSSMGroup<event_t>::SelectEvent()
{
	double u = 1.-Random::rand_uniform();
	event_t *event;
	// find the biggest event for that u < \sum omega_i
	for(typename std::vector<event_t*>::iterator e=_events.begin();
			e!=_events.end(); ++e) {
			u-=(*e)->Rate();
			if(u<=0) return *e;
	}
	// should never happen
	return _events.back();
}

}}

#endif
