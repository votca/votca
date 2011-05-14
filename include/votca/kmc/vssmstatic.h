#ifndef __VOTCA_KMC_VSSMSTATIC_H_
#define __VOTCA_KMC_VSSMSTATIC_H_

#include <vector>
#include "vssmgroup.h"

namespace votca { namespace kmc {


/**
 * \brief a specialized implementation of VSSM group for non-changing rates
 *
 * If rates do not change for all events in the group. The accumulated rate can be precalculated which allows a binary search.
 * Using this class should be used instead of VSSMGroup wherever possible.
 */
template<typename event_t>
class VSSMStatic : public VSSMGroup<event_t> {
public:
	VSSMStatic() {}

	double Rate() { return _acc_rate.back(); };

protected:

	void onEventAdded(event_t *event);

	event_t *SelectEvent();

	// precalculated accumulated rate, this allows for quick binary search
	std::vector<double> _acc_rate;
};

template<typename event_t>
inline void VSSMStatic<event_t>::onEventAdded(event_t *event)
{
	_acc_rate.push_back(_acc_rate.back() + event->Rate());

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

	while(imin - imax > 1) {
		int imid=(int)((imin+imax)*0.5);
		if(u<_acc_rate[imid])
			imin=imid;
		else
			imax=imid;
	}

	return VSSMGroup<event_t>::getEvent(imin);
}

}}

#endif

