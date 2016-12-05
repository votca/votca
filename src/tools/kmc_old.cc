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

#include <votca/xtp/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/application.h>
#include <votca/tools/tokenizer.h>
using namespace std;
using namespace votca::xtp;

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

class KMCApplication
	: public Application
{
public:
	    virtual string ProgramName() { return "kmc_run"; }
	    virtual void HelpText(std::ostream &out) { };
	    void Initialize();
	    bool EvaluateOptions();
	    void Run(void);
protected:
	    void LoadGraph(void);
		void RunKMC(void);
		void WriteOcc(void);

		map<int , node_t *> _nodes_lookup;
		vector<node_t *> _nodes;
                map<int , node_t *> _injection_lookup;
		vector<node_t *> _injection;
                string _injection_name;
		double _runtime;
		double _dt;
};

void KMCApplication::Initialize()
{
	AddProgramOptions("KMC Options")
		("graph,g", boost::program_options::value<string>(), "file which contains graph")
		("time,t", boost::program_options::value<double>(), "time to run")
	    ("seed,s", boost::program_options::value<int>(), "seed for kmc")
		("dt", boost::program_options::value<double>()->default_value(1e-4), "output frequency during run")
          ("injection,i", boost::program_options::value<string>()->default_value("*"), "name pattern for conjugated segments as injection points");
}

bool KMCApplication::EvaluateOptions()
{
	CheckRequired("graph");
	CheckRequired("seed");
	CheckRequired("time");
	_runtime = OptionsMap()["time"].as<double>();
	_dt = OptionsMap()["dt"].as<double>();
        _injection_name=OptionsMap()["injection"].as<string>();
        //cout <<_injection_name<<endl;
        //_injection_name=OptionsMap()["i"].as<string>();
	return true;
}

void KMCApplication::Run(void)
{
    srand(OptionsMap()["seed"].as<int>());
    Random::init(rand(), rand(), rand(), rand());
    LoadGraph();
    RunKMC();
}

void KMCApplication::LoadGraph() {

    Database db;
    db.Open(OptionsMap()["graph"].as<string > ());

    Statement *stmt = db.Prepare("SELECT _id,name FROM conjsegs;");

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
    cout << "Nodes: " << _nodes.size() << endl;
    cout << "Injection Points: " << _injection.size() << endl;

    delete stmt;


    int links = 0;
    stmt = db.Prepare("SELECT conjseg1, conjseg2, rate12, rate21, r_x, r_y, r_z FROM pairs;");
    while (stmt->Step() != SQLITE_DONE) {
        node_t *n1 = _nodes_lookup[stmt->Column<int>(0)];
        node_t *n2 = _nodes_lookup[stmt->Column<int>(1)];
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        vec r = vec(stmt->Column<double>(4), stmt->Column<double>(5), stmt->Column<double>(6));
        n1->AddEvent(new link_t(n2, rate12, r));
        n2->AddEvent(new link_t(n1, rate21, -r));
        links += 2;
    }
    delete stmt;
    cout << "Links: " << links << endl;
}

void KMCApplication::RunKMC(void)
{
	double t = 0;
	current=_injection[Random::rand_uniform_int(_injection.size())];
        cout <<"starting simulation at node: "<<current->_id-1<<endl;
	double next_output = _dt;
    int i=0;
    while(t<_runtime) {
    	t+=current->WaitingTime();
		current->onExecute();
    	if(t > next_output) {
    		next_output = t + _dt;
    		cout << t << ": " << r << endl;
    	}
    }
    _runtime = t;
    WriteOcc();
    cout << std::scientific << "\nKMC run finished\nAverage velocity (m/s): " << r/t*1e-9 << endl;
}

void KMCApplication::WriteOcc()
{
    Database db;
	db.Open(OptionsMap()["graph"].as<string>());
	db.Exec("BEGIN;");
	Statement *stmt = db.Prepare("UPDATE conjsegs SET occ = ? WHERE _id = ?;");
	for(int i=0; i<_nodes.size(); ++i) {
		stmt->Reset();
		stmt->Bind(1, _nodes[i]->_occ/_runtime);
		stmt->Bind(2, _nodes[i]->_id);
		stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}

int main(int argc, char **argv)
{
	KMCApplication app;
	return app.Exec(argc, argv);
}
 
