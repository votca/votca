#include <votca/kmc/vssmstatic.h>
#include <vector>
#include <iostream>


using namespace std;
using namespace votca::kmc;

struct link_t;

class node_t : public VSSMStatic<link_t> {
	double _occ;
};

struct link_t {
	double _rate;
	node_t *_dest;
};

int main(int argc, char **argv)
{
    vector<node_t *> nodes;
    return 0;
}
 
