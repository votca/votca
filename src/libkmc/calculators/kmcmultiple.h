/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *vector<Node*> 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_KMC_MULTIPLE_H
#define	__VOTCA_KMC_MULTIPLE_H

// #include <votca/kmc/vssmgroup.h>
#include <vector>
// #include <map>
#include <iostream>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <votca/tools/vec.h>
#include <votca/tools/statement.h>
#include <votca/tools/database.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include <votca/tools/random2.h>
// #include "node.h"

using namespace std;


namespace votca { namespace kmc {

static int verbose = 0; // 0=minimal output, 1=verbose output

struct Event
{
    // Event() : dr(3) {}
    int destination;
    double rate;
    votca::tools::vec::vec dr;
};


class Node
{
    public:
        int id;
        int occupied;
        int injectable;
        double escaperate;
        double occupationtime;
        vector<Event> event;
    
        double EscapeRate();
        // double CurrentEscapeRate(Node node2);
        // void Jump(Node node2);
        void AddEvent(int seg2, double rate12, votca::tools::vec::vec dr);
        void InitEscapeRate();
};


void Node::AddEvent(int seg2, double rate12, votca::tools::vec::vec dr)
{
    Event newEvent;
    newEvent.destination = seg2;
    newEvent.rate = rate12;
    newEvent.dr = dr;
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

class Chargecarrier
{
    public:
        int position;
        int id;
        Node *node;
        votca::tools::vec::vec dr_travelled;
};



void progressbar(double fraction)
{
    int totalbars = 50;
    std::cout << "\r";
    for(double bars=0; bars<double(totalbars); bars++)
    {
        if(bars<=fraction*double(totalbars))
        {
            std::cout << "|";
        }
        else
        {
            std::cout << "-";
        }
    }
    std::cout << "  " << int(fraction*1000)/10. <<" %   ";
    std::cout << std::flush;
    if(fraction*100 == 100)
    {
        std::cout << std::endl;
    }
}

int OMPinfo() 
{
    int nthreads=1, tid=0, procs, inpar=0;
    printf("\n||\n|| openMP PARALLEL COMPUTATION STATUS:\n");
    #pragma omp parallel private(tid)
    {
        # ifdef _OPENMP
        tid = omp_get_thread_num();
        if (tid == 0) 
        {
            procs = omp_get_num_procs();
            nthreads = omp_get_num_threads();
            inpar = omp_in_parallel();
            printf("|| Number of processors = %d\n", procs);
            printf("|| Number of threads = %d\n", nthreads);
    
        }
        # endif
        if (tid == 0) 
        {
            printf("|| In parallel? = %d\n||\n", inpar);
        }
}
    return nthreads;
}


class KMCMultiple : public KMCCalculator 
{
public:
    KMCMultiple() {};
   ~KMCMultiple() {};

    void Initialize(const char *filename, Property *options );
    bool EvaluateFrame();

protected:
	    vector<Node*>  LoadGraph();
            vector<double> RunVSSM(vector<Node*> node, double runtime, unsigned int numberofcharges, votca::tools::Random2 *RandomVariable);
            void WriteOcc(vector<double> occP, vector<Node*> node);

            // map<int , node_t *> _nodes_lookup;
            // vector<node_t *> _nodes;
            // map<int , node_t *> _injection_lookup;
            // vector<node_t *> _injection;
            string _injection_name;
            double _runtime;
            // double _dt;
            int _seed;
            int _numberofcharges;
            int _allowparallel;
            string _filename; // HACK
            
};

void KMCMultiple::Initialize(const char *filename, Property *options )
{
    	if (options->exists("options.kmcmultiple.runtime")) {
	    _runtime = options->get("options.kmcmultiple.runtime").as<double>();
	}
	else {
	    throw runtime_error("Error in kmcmultiple: total run time is not provided");
        }

//    	if (options->exists("options.kmcmultiple.outputtime")) {
//	    _dt = options->get("options.kmcmultiple.outputtime").as<double>();
//	}
//	else {
//	    throw runtime_error("Error in kmcmultiple: output frequency is not provided");
//        }

    	if (options->exists("options.kmcmultiple.seed")) {
	    _seed = options->get("options.kmcmultiple.seed").as<int>();
	}
	else {
	    throw runtime_error("Error in kmcmultiple: seed is not provided");
        }
        
    	if (options->exists("options.kmcmultiple.numberofcharges")) {
	    _numberofcharges = options->get("options.kmcmultiple.numberofcharges").as<int>();
	}
	else {
	    throw runtime_error("Error in kmcmultiple: number of charges is not provided");
        }

        if (options->exists("options.kmcmultiple.injection")) {
	    _injection_name = options->get("options.kmcmultiple.injection").as<string>();
	}
        else {
	    throw runtime_error("Error in kmcmultiple: injection pattern is not provided");
        }

        if (options->exists("options.kmcmultiple.allowparallel")) {
	    _allowparallel = options->get("options.kmcmultiple.allowparallel").as<int>();
	}
        else {
	    cout << "WARNING in kmcmultiple: You did not specify if parallel computation is allowed. It will be disabled." << endl;
            _allowparallel = 0;
        }

        _filename = filename;

       //cout << _seed << endl;
       //srand(_seed);
       //votca::tools::Random::init(rand(), rand(), rand(), rand());

}

vector<Node*> KMCMultiple::LoadGraph()
{
    string injectionpattern = "*";
    
    vector<Node*> node;
    
    // Load nodes
    votca::tools::Database db;
    db.Open( _filename );
    if(verbose == 1) {cout << "LOADING GRAPH" << endl << "database file: " << _filename << endl; }
    votca::tools::Statement *stmt = db.Prepare("SELECT id-1, name FROM segments;");

    int i=0;
    while (stmt->Step() != SQLITE_DONE)
    {
        Node *newNode = new Node();
        node.push_back(newNode);

        int newid = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        node[i]->id = newid;
        if (votca::tools::wildcmp(injectionpattern.c_str(), name.c_str()))
        {
            node[i]->injectable = 1;
        }
        else
        {
            node[i]->injectable = 0;
        }
        i++;
    }
    delete stmt;
    if(verbose == 1) { cout << "segments: " << node.size() << endl; }
    
    // Load pairs and rates
    int numberofpairs = 0;
    stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', rate12e AS 'rate', drX, drY, drZ FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', rate21e AS 'rate', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ' FROM pairs ORDER BY segment1;");
    while (stmt->Step() != SQLITE_DONE)
    {
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);
        double rate12 = stmt->Column<double>(2);
        votca::tools::vec::vec dr = votca::tools::vec::vec(stmt->Column<double>(3), stmt->Column<double>(4), stmt->Column<double>(5));
        node[seg1]->AddEvent(seg2,rate12,dr);
        numberofpairs ++;
    }    
    delete stmt;

    if(verbose == 1) { cout << "pairs: " << numberofpairs/2 << endl; }
    
    // Calculate initial escape rates
    for(unsigned int i=0; i<node.size(); i++)
    {
        node[i]->InitEscapeRate();
    }
    return node;
}

vector<double> KMCMultiple::RunVSSM(vector<Node*> node, double runtime, unsigned int numberofcharges, votca::tools::Random2 *RandomVariable)
{
    cout << endl << "Algorithm: VSSM for Multiple Charges" << endl;
    cout << "number of charges: " << numberofcharges << endl;
    cout << "number of nodes: " << node.size() << endl;
    
    if(numberofcharges > node.size())
    {
        throw runtime_error("ERROR in kmcmultiple: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
    }
    
    double outputfrequency = runtime/100;
    vector<double> occP(node.size(),0.);

    if(numberofcharges>node.size())
    {
        throw runtime_error("Error in kmcstandalone: Your number of charges is larger than the number of nodes. This conflicts with single occupation.");
    }

    // Injection
    vector< Chargecarrier* > carrier;
    for (unsigned int i = 0; i < numberofcharges; i++)
    {
        Chargecarrier *newCharge = new Chargecarrier;      
        newCharge->id = i;
        newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];
        while(newCharge->node->occupied == 1 || newCharge->node->injectable != 1)
        {   // maybe already occupied? or maybe not injectable?
            newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];
        }
        newCharge->node->occupied = 1;
        newCharge->dr_travelled = votca::tools::vec::vec(0,0,0);
        cout << "starting position for charge " << i+1 << ": segment " << newCharge->node->id+1 << endl;
        carrier.push_back(newCharge);
    }    

    double time = 0.;
    int step = 0;
    double nextoutput = outputfrequency;
    
    progressbar(0.);
    while(time < runtime)
    {
        
        // determine which electron will escape
        Node* do_oldnode;
        Node* do_newnode;
        Chargecarrier* do_affectedcarrier;
        double maxprob = 0.;
        double newprob = 0.;
        votca::tools::vec::vec dr;
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            newprob = carrier[i]->node->EscapeRate() * (1-RandomVariable->rand_uniform());
            if(newprob > maxprob)
            {
                maxprob = newprob;
                do_oldnode = carrier[i]->node;
                do_affectedcarrier = carrier[i];
            }
        }
        if(verbose >= 1) {cout << "Charge number " << do_affectedcarrier->id+1 << " which is sitting on segment " << do_oldnode->id+1 << " will escape!" << endl ;}
        
        // determine where it will jump to
        if(verbose >= 1) {cout << "There are " << do_oldnode->event.size() << " possible jumps for this charge:"; }
        maxprob = 0;
        for(unsigned int j=0; j<do_oldnode->event.size(); j++)
        {
            if(verbose >= 1) { cout << " " << do_oldnode->event[j].destination ; }
            newprob = do_oldnode->event[j].rate * (1-RandomVariable->rand_uniform());
            if(newprob > maxprob)
            {
                maxprob = newprob;
                do_newnode = node[do_oldnode->event[j].destination];
                dr = do_oldnode->event[j].dr;
            }
        }
        
        if(verbose >= 1) {cout << "." << endl;}

        // go forward in time
        double accumulatedrate = 0;
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            accumulatedrate += carrier[i]->node->EscapeRate();
        }
        double dt = 0;
        if(accumulatedrate == 0)
        {   // this should not happen: no possible jumps defined for a node
            throw runtime_error("Error in kmcsingle: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
        }
        else
        {

            double rand_u = 1-RandomVariable->rand_uniform();
            while(rand_u == 0)
            {
                cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                rand_u = 1-RandomVariable->rand_uniform();
            }

                
            dt = -1 / accumulatedrate * log(rand_u);
        }
        time += dt;
        step += 1;
        
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            carrier[i]->node->occupationtime += dt;
        }
        
        
        // if the new segment is unoccupied: jump; if not: nothing??
        if(do_newnode->occupied == 0)
        {
            do_newnode->occupied = 1;
            do_oldnode->occupied = 0;
            do_affectedcarrier->node = do_newnode;
            do_affectedcarrier->dr_travelled += dr;
            if(verbose >= 1) {cout << "Charge has jumped to segment: " << do_newnode->id+1 << "." << endl << endl;}
        }
        else
        {
            if(verbose >= 1) {cout << "Selected segment: " << do_newnode->id+1 << " is alread occupied. No jump." << endl << endl;}
        }
        
        if(time > nextoutput)
        {
            nextoutput = time + outputfrequency;
            progressbar(time/runtime);
        }
    }
    progressbar(1.);

    
    // calculate occupation probabilites from occupation times    
    for(unsigned int j=0; j<node.size(); j++)
    {   
        occP[j] = node[j]->occupationtime / time;
    }

    cout << endl << "finished KMC simulation after " << step << " steps." << endl;
    votca::tools::vec::vec dr_travelled = votca::tools::vec::vec (0,0,0);
    cout << endl << "Average velocities: " << endl;
    for(unsigned int i=0; i<numberofcharges; i++)
    {
        cout << std::scientific << "    charge " << i+1 << ": " << carrier[i]->dr_travelled/time*1e-9 << endl;
        dr_travelled += carrier[i]->dr_travelled;
    }
    dr_travelled /= numberofcharges;
    cout << std::scientific << "  Overall average velocity (m/s): " << dr_travelled/time*1e-9 << endl << endl;
    return occP;
}


void KMCMultiple::WriteOcc(vector<double> occP, vector<Node*> node)
{
    votca::tools::Database db;
    cout << "Opening for writing " << _filename << endl;
	db.Open(_filename);
	db.Exec("BEGIN;");
	votca::tools::Statement *stmt = db.Prepare("UPDATE segments SET occPe = ? WHERE id = ?;");  // electron occ. prob., check (think about) this
	for(unsigned int i=0; i<node.size(); ++i)
        {
	    stmt->Reset();
	    stmt->Bind(1, occP[i]);
	    stmt->Bind(2, node[i]->id+1);
	    stmt->Step();
	}
	db.Exec("END;");
	delete stmt;
}

bool KMCMultiple::EvaluateFrame()
{
    // double runtime = 1E6;
    // int seed  = 23;
    // unsigned int numberofcharges = 2;
    
    unsigned int numberofthreads = 1;
    if(_allowparallel == 1)
    {
        numberofthreads = OMPinfo();
    }
    

    std::cout << "-----------------------------------" << std::endl;      
    std::cout << "      KMC FOR MULTIPLE CHARGES" << std::endl;
    std::cout << "-----------------------------------" << std::endl << std::endl;      
 
    
    
    // Initialise random number generator
    // each thread i in a parallel computation needs is own set RandomVariable[i]
    if(verbose == 1) { cout << endl << "Initialising random number generator" << endl; }
    srand(_seed); // srand expects any integer in order to initialise the random number generator
    vector< votca::tools::Random2* > RandomVariable;
    for (unsigned int i = 0; i < numberofthreads; i++)
    {
        votca::tools::Random2 *newRandomVariable = new votca::tools::Random2();      
        RandomVariable.push_back(newRandomVariable);
        RandomVariable[i]->init(rand(), rand(), rand(), rand());
    }
    
    // VSSM KMC algorithm
    //cout << endl << "KMC SIMULATION" << endl;
    cout << "number of threads: " << numberofthreads << endl;

    vector<Node*> node;
    node = KMCMultiple::LoadGraph();
    vector<double> occP(node.size(),0.);
    vector< vector< double > > occPOneRun ( numberofthreads, vector<double> ( node.size(), 0. ) );

    # ifdef _OPENMP
    (void) omp_set_num_threads(numberofthreads);
    # endif
    #pragma omp parallel private(node)
    {
        node = KMCMultiple::LoadGraph();
        int thread = 0;
        # ifdef _OPENMP
        thread = omp_get_thread_num();
        # endif
        occPOneRun[thread] = KMCMultiple::RunVSSM(node, _runtime/numberofthreads, _numberofcharges, RandomVariable[thread]);
    }

    
    // get mean of multiple runs
    for(unsigned int j=0; j<occP.size();j++) 
    {
        for(unsigned int thread=0; thread<numberofthreads; thread++)
        {
            occP[j] += occPOneRun[thread][j];
        }
        occP[j] /= numberofthreads;
    }
    
    // output occupation probabilites
    for(unsigned int thread=0; thread<numberofthreads; thread++)
    {
        for(unsigned int j=0; j<occPOneRun[thread].size();j++) 
        {
            if(occPOneRun[thread][j] > 0)
            {
                cout << "[thread " << thread+1 << "] "<<"occupation probability " << node[j]->id+1 << ": " << occPOneRun[thread][j] << endl;
            }
        }
    }
    // output occupation probabilites
    for(unsigned int j=0; j<occP.size();j++) 
    {
        if(occP[j] > 0)
        {
            cout << "occupation probability " << node[j]->id+1 << ": " << occP[j] << endl;
        }
    }
    
    // write occupation probabilites
    KMCMultiple::WriteOcc(occP, node);
    
    return true;
}


}}


#endif	/* __VOTCA_KMC_MULTIPLE_H */
