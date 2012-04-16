/* 
 * File:   node.h
 * Author: malafeev
 *
 * Created on November 29, 2011, 3:17 PM
 */

#ifndef NODE_H
#define	NODE_H

using namespace std;
using namespace votca::kmc;

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

#endif	/* NODE_H */

