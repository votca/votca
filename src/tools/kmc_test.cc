#include <votca/kmc/vssmgroup.h>
#include <vector>
#include <iostream>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::kmc;

struct link_t;

class node_t : public VSSMGroup<link_t> {
	double _occ;
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
    vector<node_t *> nodes;
    Random::init(1, 2, 3, 4);
    for(int x=0; x<10; ++x)
        for(int y=0; y<10; ++y)
            for(int z=0; z<10; ++z)
            	nodes.push_back(new node_t());
    for(int x=0; x<10; ++x)
        for(int y=0; y<10; ++y)
            for(int z=0; z<10; ++z)
				for(int dx=-1; dx<=1; ++dx)
					for(int dy=-1; dy<=1; ++dy)
						for(int dz=-1; dz<=1; ++dz) {
							if(dx==0 && dy==0 && dz==0)
								continue;
							node_t *n = nodes[x*100 + y*10 + z];
							node_t *dest = nodes[((x+dx+10)%10)*100 + ((y+dy+10)%10)*10 + (z+dz+10)%10];
							n->AddEvent(new link_t(dest, 1.0, vec(dx, dy, dz)));
						}
    current=nodes[0];
    for(int n=0; n<100; ++n) {
    	for(int m=0; m<100000; ++m) {
    		current->onExecute();
    	}
    	cout << n << " " << r << endl;
    }
    return 0;
}
 
