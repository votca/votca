#include <votca/kmc/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>

using namespace std;
using namespace votca::kmc;

struct link_t;

class node_t : public VSSMGroup<link_t> {
  public:
	node_t(int id)
	  : _id(id) {}
	double _occ;
	int _id;
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

int main(int argc, char **argv)
{
	map<int , node_t *> nodes;
    Database db;
	db.Open("state.db");
	
    Random::init(1, 2, 3, 4);

	Statement *stmt = db.Prepare("SELECT id FROM crgunits;");
	
	while(stmt->Step() != SQLITE_DONE) {
		int id = stmt->Column<int>(0);
		nodes[id] = new node_t(id);
	}
	delete stmt;
	cout << "Nodes: " << nodes.size() << endl;
	
	int links=0;
	stmt = db.Prepare("SELECT crgunit1, crgunit2, rate12, rate21, r_x, r_y, r_z FROM pairs;");
	while(stmt->Step() != SQLITE_DONE) {
	  node_t *n1 = nodes[stmt->Column<int>(0)];
	  node_t *n2 = nodes[stmt->Column<int>(1)];
	  double rate12 = stmt->Column<double>(2);
	  double rate21 = stmt->Column<double>(3);
	  vec r = vec(stmt->Column<double>(4),stmt->Column<double>(5),stmt->Column<double>(6));
	  n1->AddEvent(new link_t(n2, rate12, r));
	  n2->AddEvent(new link_t(n1, rate21, -r));
	  links+=2;
	}
	delete stmt;
	cout << "Links: " << links << endl;

	double t = 0;
	current=nodes[10];
    while(t<1e-1) {
	  for(int m=0; m<2000000; ++m) {
		    t+=current->WaitingTime();
    		current->onExecute();
    	}
    	cout << t << ": " << r << endl;
    }
    return 0;
}
 
