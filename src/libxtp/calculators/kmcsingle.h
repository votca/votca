/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_SINGLE_H
#define	__VOTCA_KMC_SINGLE_H

#include <votca/xtp/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>

#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include <votca/tools/random.h>

namespace votca { namespace xtp {
    
using namespace std;
using namespace votca::xtp;
using namespace votca::tools;

struct link_t;
struct event_t;

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
	double WaitingTime() { return _waiting_time; }

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
	//double max = Rate();
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



class node_t : public VSSMGroup<link_t> {
  public:
	node_t(int id)
	  : _occ(0),_id(id) {}
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
	link_t(double rate, node_t *dest, vec r)
	:  _rate(rate),_dest(dest), _r(r) {}
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

class KMCSingle : public KMCCalculator 
{
public:
    KMCSingle() {};
   ~KMCSingle() {};

    void Initialize(const char *filename, Property *options , const char *outputfile);
    bool EvaluateFrame();

protected:
	    void LoadGraph();
            void RunKMC(void);
            void WriteOcc(void);

            map<int , node_t *> _nodes_lookup;
            vector<node_t *> _nodes;
            map<int , node_t *> _injection_lookup;
            vector<node_t *> _injection;
            string _injection_name;
            double _runtime;
            double _dt;
            int _seed;
            string _filename; // HACK
            string _outputfile;
            
};

void KMCSingle::Initialize(const char *filename, Property *options , const char *outputfile)
{
    	if (options->exists("options.kmcsingle.runtime")) {
	    _runtime = options->get("options.kmcsingle.runtime").as<double>();
	}
	else {
	    throw invalid_argument("Error in kmcsingle: total run time is not provided");
        }

    	if (options->exists("options.kmcsingle.outputtime")) {
	    _dt = options->get("options.kmcsingle.outputtime").as<double>();
	}
	else {
	    throw invalid_argument("Error in kmcsingle: output frequency is not provided");
        }

    	if (options->exists("options.kmcsingle.seed")) {
	    _seed = options->get("options.kmcsingle.seed").as<int>();
	}
	else {
	    throw invalid_argument("Error in kmcsingle: seed is not provided");
        }
        
   	if (options->exists("options.kmcsingle.injection")) {
	    _injection_name = options->get("options.kmcsingle.injection").as<string>();
	}
        else {
	    throw invalid_argument("Error in kmcsingle: injection pattern is not provided");
        }
        
        _filename = filename;
        _outputfile = outputfile;

       //cout << _seed << endl;
       srand(_seed);
       votca::tools::Random::init(rand(), rand(), rand(), rand());

}

bool KMCSingle::EvaluateFrame()
{

    LoadGraph();
    RunKMC();
    return true;
}

void KMCSingle::LoadGraph() {

    Database db;
    db.Open( _filename );
    cout << " Loading graph from " << _filename << endl;
    Statement *stmt = db.Prepare("SELECT _id, name FROM segments;");

    while (stmt->Step() != SQLITE_DONE) {
        int id = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        node_t *n =new node_t(id);
        _nodes.push_back(n);
        _nodes_lookup[id] = _nodes.back();
        if (wildcmp(_injection_name.c_str(), name.c_str())) {
            _injection.push_back(n);
            _injection_lookup[id] = _injection.back();
        }
    }
    //delete stmt;
    cout << "  -Nodes: " << _nodes.size() << endl;
    cout << "seed:" << _seed << endl;
    cout << "  -Injection Points: " << _injection.size() << endl;

    delete stmt;

    int links = 0;
    stmt = db.Prepare("SELECT seg1, seg2, rate12e, rate21e, drX, drY, drZ FROM pairs;"); // electron rates, kmcsingle still lacks capability to treat holes
    while (stmt->Step() != SQLITE_DONE) {
        node_t *n1 = _nodes_lookup[stmt->Column<int>(0)];
        node_t *n2 = _nodes_lookup[stmt->Column<int>(1)];
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        vec r = vec(stmt->Column<double>(4), stmt->Column<double>(5), stmt->Column<double>(6));
        n1->AddEvent(new link_t(rate12,n2, r));
        n2->AddEvent(new link_t( rate21,n1, -r));
        links += 2;
        
        if ( votca::tools::globals::verbose ) {
            cout << "rate12=" << rate12 << endl;
            cout << "rate21=" << rate21 << endl;
            cout << "r=" << r << endl;
        }
    }
    delete stmt;
    cout << "  -Links: " << links << endl;

}

void KMCSingle::RunKMC(void)
{

	double t = 0;

        srand(_seed);
        Random::init(rand(), rand(), rand(), rand());

        // cout << " seed:size:site " << _seed << ":" << _injection.size() << ":" << Random::rand_uniform_int(_injection.size()) << endl;
	current=_injection[Random::rand_uniform_int(_injection.size())];
        cout <<" Starting simulation at node: "<<current->_id-1<<endl;
        if (current->Rate() == 0) {
            int id = current->_id;
            throw std::runtime_error("Injected on an unconnected site, id = "+boost::lexical_cast<string>(id));            
        }
	double next_output = _dt;
    //int i=0;
    while(t<_runtime) {
    	t+=current->WaitingTime();
		current->onExecute(); // this line causes a Segmentation fault
    	if(t > next_output) {
    		next_output = t + _dt;
    		cout << t << ": " << r << endl;
    	}
    }
    _runtime = t;
    WriteOcc();
    cout << std::scientific << "\nKMC run finished\nAverage velocity (m/s): " << r/t*1e-9 << endl;
}

void KMCSingle::WriteOcc()
{
    Database db;
    cout << "Opening for writing " << _filename << endl;
	db.Open(_filename);
	db.Exec("BEGIN;");
	Statement *stmt = db.Prepare("UPDATE segments SET occPe = ? WHERE _id = ?;");  // electron occ. prob., check (think about) this
	for(unsigned i=0; i<_nodes.size(); ++i) {
		stmt->Reset();
		stmt->Bind(1, _nodes[i]->_occ/_runtime);
		stmt->Bind(2, _nodes[i]->_id);
		stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}

}}
#endif	/* __VOTCA_KMC_SINGLE_H */
