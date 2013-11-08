/* 
 * author: Kordt
 */

#ifndef NODE_H
#define	NODE_H
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::kmc;


// KMCMULTIPLE PART //
typedef votca::tools::vec myvec;

struct Event
{
    int destination;
    double rate;
    votca::tools::vec dr;
    
    // stuff for Coulomb interaction
    double Jeff2;
    double reorg_out;
    double initialrate;
};


class Node
{
    public:
        int id;
        int occupied;
        int injectable;
        double escaperate;
        double occupationtime;
        myvec position;
        vector<Event> event;
        // stuff for Coulomb interaction:
        double siteenergy;
        double reorg_intorig; // UnCnN
        double reorg_intdest; // UcNcC
    
        double EscapeRate();
        void AddEvent(int seg2, double rate12, myvec dr, double Jeff2, double reorg_out);
        void InitEscapeRate();
};


void Node::AddEvent(int seg2, double rate12, myvec dr, double Jeff2, double reorg_out)
{
    Event newEvent;
    newEvent.destination = seg2;
    newEvent.rate = rate12;
    newEvent.initialrate = rate12;
    newEvent.dr = dr;
    newEvent.Jeff2 = Jeff2;
    newEvent.reorg_out = reorg_out;
    this->event.push_back(newEvent);
}


void Node::InitEscapeRate()
{
    double newEscapeRate = 0;
    for(unsigned int i=0; i<this->event.size();i++)
    {
        newEscapeRate += this->event[i].rate;
    }
    this->escaperate = newEscapeRate;
    // cout << "Escape rate for segment " << this->id << " was set to " << newEscapeRate << endl;
}


double Node::EscapeRate()
{
    return escaperate;
}
// END KMCMULTIPLE PART //




// KMCSINGLE PART //
/*
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

#endif	/* NODE_H */

