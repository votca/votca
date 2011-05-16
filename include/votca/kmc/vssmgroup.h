#ifndef __VOTCA_KMC_VSSMGROUP_H_
#define __VOTCA_KMC_VSSMGROUP_H_

#include <votca/tools/random.h>
#include <vector>
#include <iostream>

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
 * \todo add a refresh function if rates have changed
 */
template<typename event_t>
class VSSMGroup {
public:
	VSSMGroup() {
		_acc_rate.push_back(0);
	}

	/***
	 * \brief add event to list of possible events
	 */
	void AddEvent(event_t *event) {
		_events.push_back(event);
		_acc_rate.push_back(_acc_rate.back() + event->Rate());
		UpdateWaitingTime();
	}

	/**
	 * \brief get the total rate of this group
	 */
	double Rate() { return _acc_rate.back(); };

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
	 * \brief select an event to execute based on linear search O(N)
	 */
	event_t *SelectEvent_LinearSearch();

	/**
	 * \brief select an event to execute based on binary search O(log N)
	 */
	event_t *SelectEvent_BinarySearch();

	std::vector<event_t*> _events;
	std::vector<double> _acc_rate;
	double _waiting_time;
};

template<typename event_t>
inline void VSSMGroup<event_t>::onExecute()
{
	SelectEvent_BinarySearch()->onExecute();
	UpdateWaitingTime();
}

template<typename event_t>
inline event_t *VSSMGroup<event_t>::SelectEvent_LinearSearch()
{
	double u = (1.-Random::rand_uniform())*Rate();
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

template<typename event_t>
event_t *VSSMGroup<event_t>::SelectEvent_BinarySearch()
{
	double u = 1.-Random::rand_uniform();
	u=u*Rate();
	double max = Rate();
	// to a binary search in accumulated events
	int imin=0;
	int imax=_acc_rate.size();

	while(imax - imin > 1) {
		int imid=(int)((imin+imax)*0.5);
		//std::cout << u << " " << _acc_rate[imid] << std::endl;
		if(u<=_acc_rate[imid])
			imax=imid;
		else
			imin=imid;
	}
   // std::cout << imin << " " << Rate() <<" " << u << std::endl;
	return _events[imin];
}

}}

#endif
