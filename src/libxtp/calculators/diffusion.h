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

#ifndef __VOTCA_DIFFUSION_H
#define	__VOTCA_DIFFUSION_H

#include <votca/xtp/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/matrix.h>
#include "node.h"

using namespace std;
using namespace votca::xtp;
using namespace votca::tools;



class Diffusion : public KMCCalculator
{
public:
    Diffusion() {zero_border=0.0;};
   ~Diffusion() {};

    void Initialize(const char *filename, Property *options , const char *outputfile);
    bool EvaluateFrame();
    double zero_border;
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
            matrix diffusion; //matrix for diffusion tensor
            int    number_of_points; // number of vectors in diffusion tensor calculation
            matrix::eigensystem_t diff_tensor_eigensystem;
            
};

void Diffusion::Initialize(const char *filename, Property *options , const char *outputfile)
{
    	if (options->exists("options.diffusion.runtime")) {
	    _runtime = options->get("options.diffusion.runtime").as<double>();
	}
	else {
	    std::runtime_error("Error in diffusion: total run time is not provided");
        }

    	if (options->exists("options.diffusion.outputtime")) {
	    _dt = options->get("options.diffusion.outputtime").as<double>();
	}
	else {
	    std::runtime_error("Error in diffusion: output frequency is not provided");
        }

    	if (options->exists("options.diffusion.seed")) {
	    _seed = options->get("options.diffusion.seed").as<int>();
	}
	else {
	    std::runtime_error("Error in diffusion: seed is not provided");
        }
        
   	if (options->exists("options.diffusion.injection")) {
	    _injection_name = options->get("options.diffusion.injection").as<string>();
	}
	else {
	    std::runtime_error("Error in diffusion: injection pattern is not provided");
        }    

        if(options->exists("options.rates.field"))
        {
            
            bool  field_is_zero;
            field_is_zero = true;

            if( options->get("options.rates.field").as<vec>().x()>= zero_border )
                {
                field_is_zero = false;
                }
            if( options->get("options.rates.field").as<vec>().y()>= zero_border )
                {
                field_is_zero = false;
                }
            if( options->get("options.rates.field").as<vec>().z()>= zero_border )
                {
                field_is_zero = false;
                }

            if(!field_is_zero)
                {
                cout<<"WARNING: Electric field is not zero!"<<endl;
                }
        }else
        {
               cout<<"WARNING: Can't find electric field in the option file"<<endl;
               cout<<"check whether it is zero"<<endl;
        }

        
        _filename = filename;
        _outputfile = outputfile;
}

bool Diffusion::EvaluateFrame()
{
    srand(_seed);
    Random::init(rand(), rand(), rand(), rand());
    LoadGraph();
    RunKMC();
    return true;
}

void Diffusion::LoadGraph() {

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
    cout << "  -Injection Points: " << _injection.size() << endl;

    delete stmt;


    int links = 0;
    stmt = db.Prepare("SELECT seg1, seg2, rate12e, rate21e, drX, drY, drZ FROM pairs;");
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
    cout << "  -Links: " << links << endl;
}

void Diffusion::RunKMC(void)
{
    double t = 0;
        number_of_points = 0;
        //put zeros in diffusion tensor matrix
        diffusion.ZeroMatrix();

        srand(_seed);
        Random::init(rand(), rand(), rand(), rand());


	current=_injection[Random::rand_uniform_int(_injection.size())];
        cout <<" Starting simulation at node: "<<current->_id-1<<endl;
	double next_output = _dt;
    int i=0;
    while(t<_runtime) {
    	t+=current->WaitingTime();
		current->onExecute();
    	if(t > next_output) {
    		next_output = t + _dt;
                number_of_points++;

                
                // calculation of diffusion tensor
                diffusion += r|r;

    		cout << t << ": " << r << endl;
    	}
    }
    _runtime = t;
    WriteOcc();
    cout << std::scientific << "\nKMC run finished\nAverage velocity (m/s): " << r/t*1e-9 << endl;
    
    // calculating diffusion and printing matrix
    diffusion /= number_of_points*2*t;
    diffusion *= 1e-18;
    cout<<endl<<"Diffusion tensor:"<<endl<<diffusion<<endl;


    //solve eigensystem
    diffusion.SolveEigensystem(diff_tensor_eigensystem);
    //cout<<diff_tensor_eigensystem<<endl;


    cout<<endl<<"Eigenvalues: "<<endl<<endl;
    for(int i=0; i<=2; i++)
    {
        cout<<"Eigenvalue: "<<diff_tensor_eigensystem.eigenvalues[i]<<endl<<"Eigenvector: ";
               
        cout<<diff_tensor_eigensystem.eigenvecs[i].x()<<"   ";
        cout<<diff_tensor_eigensystem.eigenvecs[i].y()<<"   ";
        cout<<diff_tensor_eigensystem.eigenvecs[i].z()<<endl<<endl;
    }
    cout<<endl<<endl<<endl;




}

void Diffusion::WriteOcc()
{
    Database db;
    cout << "Opening for writing " << _filename << endl;
	db.Open(_filename);
	db.Exec("BEGIN;");
	Statement *stmt = db.Prepare("UPDATE segments SET occPe = ? WHERE _id = ?;");
	for(int i=0; i<_nodes.size(); ++i) {
		stmt->Reset();
		stmt->Bind(1, _nodes[i]->_occ/_runtime);
		stmt->Bind(2, _nodes[i]->_id);
		stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}

#endif	/* __VOTCA_DIFFUSION_H */
