/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_KMC_GNODE_H
#define	_VOTCA_KMC_GNODE_H

#include <votca/tools/vec.h>
#include <votca/kmc/link.h>

using namespace std;
using namespace votca::kmc;


// KMCMULTIPLE PART //
typedef votca::tools::vec myvec;


namespace votca { namespace kmc {

class GNode
{
    public:
        int id;
        int occupied;
        int injectable;
        double escaperate;
        double occupationtime;
        myvec position;
        vector<Link> event;
        // stuff for Coulomb interaction:
        double siteenergy;
        double reorg_intorig; // UnCnN
        double reorg_intdest; // UcNcC
    
        double EscapeRate();
        void AddEvent(int seg2, double rate12, myvec dr, double Jeff2, double reorg_out);
        void InitEscapeRate();
};


void GNode::AddEvent(int seg2, double rate12, myvec dr, double Jeff2, double reorg_out)
{
    Link newEvent;
    newEvent.destination = seg2;
    newEvent.rate = rate12;
    newEvent.initialrate = rate12;
    newEvent.dr = dr;
    newEvent.Jeff2 = Jeff2;
    newEvent.reorg_out = reorg_out;
    this->event.push_back(newEvent);
}


void GNode::InitEscapeRate()
{
    double newEscapeRate = 0;
    for(unsigned int i=0; i<this->event.size();i++)
    {
        newEscapeRate += this->event[i].rate;
    }
    this->escaperate = newEscapeRate;
    // cout << "Escape rate for segment " << this->id << " was set to " << newEscapeRate << endl;
}


double GNode::EscapeRate()
{
    return escaperate;
}
// END KMCMULTIPLE PART //


/*
// KMCSINGLE PART //
struct link_t;

class node_t : public VSSMGroup<link_t> {
  public:
	node_t(int id)
	  : _id(id), _occ(0) {}
	double _occ;
	int _id;

	void onExecute() {
                _occ+=WaitingTime();
		VSSMGroup<link_t>::onExecute();
	}
};

node_t *current;
vec r(0,0,0);

struct link_t {
	link_t(node_t *dest, double rate, vec r)
	: _dest(dest), _rate(rate), _r(r) {}
	double Rate() {
		return _rate;
	}

	void onExecute() {
		r+=_r;
		current = _dest;
	}
	double _rate;
	node_t *_dest;
	vec _r;
};
// END KMCSINGLE PART // 
*/

}}

#endif	/* _VOTCA_KMC_GNODE_H */

