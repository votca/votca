#ifndef __VOTCA_KMC_VSSMSTATIC_H_
#define __VOTCA_KMC_VSSMSTATIC_H_

#include <vector>
#include "vssmgroup.h"
#include <iostream>

namespace votca { namespace kmc {


/**
 * \brief a specialized implementation of VSSM group for non-changing rates
 *
 * If rates do not change for all events in the group. The accumulated rate can be precalculated which allows a binary search.
 * Using this class should be used instead of VSSMGroup wherever possible.
 */
template<typename event_t>
class VSSMStatic  {
public:
	VSSMStatic() {}
	void AddEvent(event_t *event) {
			_events.push_back(event);
			onEventAdded(event);
			UpdateWaitingTime();
		}
	double Rate() { return _acc_rate.back(); };
	void UpdateWaitingTime() {
		_waiting_time = -log( 1.0 - Random::rand_uniform() ) / Rate();
	}

	double WaitingTime();
	void onExecute();

protected:

	void onEventAdded(event_t *event);

	event_t *SelectEvent();

	// precalculated accumulated rate, this allows for quick binary search
	std::vector<double> _acc_rate;

	std::vector<event_t*> _events;
	double _waiting_time;
};

template<typename event_t>
inline void VSSMStatic<event_t>::onExecute()
{
	SelectEvent()->onExecute();
	UpdateWaitingTime();
}

template<typename event_t>
inline void VSSMStatic<event_t>::onEventAdded(event_t *event)
{
	if(_acc_rate.size()==0) {
		_acc_rate.push_back(0);
		_acc_rate.push_back(event->Rate());
	}
	else _acc_rate.push_back(_acc_rate.back() + event->Rate());
}

template<typename event_t>
event_t *VSSMStatic<event_t>::SelectEvent()
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

