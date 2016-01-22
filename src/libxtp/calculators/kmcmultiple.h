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
 * author: Kordt
 */

#ifndef __VOTCA_KMC_MULTIPLE_H
#define	__VOTCA_KMC_MULTIPLE_H

// #include <votca/xtp/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath> // needed for abs(double)

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/tools/statement.h>
#include <votca/tools/database.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include <votca/tools/random2.h>

#include <tr1/unordered_map>
#include <votca/xtp/gnode.h>
#include <math.h> // needed for fmod()

using namespace std;

namespace votca { namespace xtp {

    const static double kB   = 8.617332478E-5; // eV/K
    const static double hbar = 6.5821192815E-16; // eV*s
    const static double eps0 = 8.85418781762E-12/1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
    const static double epsr = 3.0; // relative material permittivity
    const static double Pi   = 3.14159265358979323846;

//int OMPinfo() 
//{
//    int nthreads=1, tid=0, procs, inpar=0;
//    printf("\n||\n|| openMP PARALLEL COMPUTATION STATUS:\n");
//    #pragma omp parallel private(tid)
//    {
//        # ifdef _OPENMP
//        tid = omp_get_thread_num();
//        if (tid == 0) 
//        {
//            procs = omp_get_num_procs();
//            nthreads = omp_get_num_threads();
//            inpar = omp_in_parallel();
//            printf("|| Number of processors = %d\n", procs);
//            printf("|| Number of threads = %d\n", nthreads);
//    
//        }
//        # endif
//        if (tid == 0) 
//        {
//            printf("|| In parallel? = %d\n||\n", inpar);
//        }
//    }
//    return nthreads;
//}


class KMCMultiple : public KMCCalculator 
{
public:
    KMCMultiple() {};
   ~KMCMultiple() {};

    void Initialize(const char *filename, Property *options, const char *outputfile );
    bool EvaluateFrame();



protected:
        
        
            class Chargecarrier
            {
                public:
                    int id;
                    GNode *node;
                    vec dr_travelled;
                    double tof_travelled;
                    
            };
            typedef std::tr1::unordered_map<unsigned long, double> CoulombMap;
            typedef CoulombMap::const_iterator CoulombIt;

	    vector<GNode*>  LoadGraph();
	    KMCMultiple::CoulombMap LoadCoulomb(int numberofnodes);
            vector<double> RunVSSM(vector<GNode*> node, double runtime, unsigned int numberofcharges, votca::tools::Random2 *RandomVariable, CoulombMap coulomb);
            void WriteOcc(vector<double> occP, vector<GNode*> node);
            void InitialRates(vector<GNode*> node);
            void RateUpdateCoulomb(vector<GNode*> &node,  vector< Chargecarrier* > &carrier, CoulombMap &coulomb);
            void InitBoxSize(vector<GNode*> node);
            void progressbar(double fraction);


            string _injection_name;
            string _injectionmethod;
            //int _injectionfree;
            double _runtime;
            double _maxrealtime;
            int _seed;
            int _numberofcharges;
            int _allowparallel;
            double _fieldX;
            double _fieldY;
            double _fieldZ;
            vec _field;
            double _outputtime;
            string _trajectoryfile;
            string _timefile;
            string _carriertype;
            int _explicitcoulomb;
            double _temperature;
            string _filename;
            string _outputfile;
            double _q;
            double _boxsizeX;
            double _boxsizeY;
            double _boxsizeZ;
            string _rates;
            int    _tof;
            double _toflength;
            string _tofdirection;
            double _coulomboffset;
};




void KMCMultiple::Initialize(const char *filename, Property *options, const char *outputfile )
{
    	if (options->exists("options.kmcmultiple.runtime")) {
	    _runtime = options->get("options.kmcmultiple.runtime").as<double>();
	}
	else {
	    throw runtime_error("Error in kmcmultiple: total run time is not provided");
        }
    	if (options->exists("options.kmcmultiple.maxrealtime")) {
	    _maxrealtime = options->get("options.kmcmultiple.maxrealtime").as<double>();
	}
        else{
            _maxrealtime = 1E10; // maximal real time in hours
        }
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
        
        if (options->exists("options.kmcmultiple.injectionmethod")) {
	    _injectionmethod = options->get("options.kmcmultiple.injectionmethod").as<string>();
	}
        else {
	    cout << "WARNING in kmcmultiple: You did not specify an injection method. It will be set to random injection." << endl;
            _injectionmethod = "random";
        }
        if (_injectionmethod != "random" && _injectionmethod != "equilibrated")
        {
	    cout << "WARNING in kmcmultiple: Unknown injection method. It will be set to random injection." << endl;
            _injectionmethod = "random";
            
        }

        if (options->exists("options.kmcmultiple.fieldX")) {
	    _fieldX = options->get("options.kmcmultiple.fieldX").as<double>();
	}
        else {
            _fieldX = 0;
        }
        if (options->exists("options.kmcmultiple.fieldY")) {
	    _fieldY = options->get("options.kmcmultiple.fieldY").as<double>();
	}
        else {
            _fieldY = 0;
        }
        if (options->exists("options.kmcmultiple.fieldZ")) {
	    _fieldZ = options->get("options.kmcmultiple.fieldZ").as<double>();
	}
        else {
            _fieldZ = 0;
        }
        _field = vec(_fieldX,_fieldY,_fieldZ);
        if (options->exists("options.kmcmultiple.outputtime")) {
	    _outputtime = options->get("options.kmcmultiple.outputtime").as<double>();
	}
        else {
            _outputtime = 0;
        }
        if (options->exists("options.kmcmultiple.trajectoryfile")) {
	    _trajectoryfile = options->get("options.kmcmultiple.trajectoryfile").as<string>();
	}
        else {
            _trajectoryfile = "trajectory.csv";
        }
        if (_trajectoryfile == "") {
            _trajectoryfile = "trajectory.csv";
        }
        if (options->exists("options.kmcmultiple.timefile")) {
	    _timefile = options->get("options.kmcmultiple.timefile").as<string>();
	}
        else {
            _timefile = "timedependence.csv";
        }
        if (_timefile == "") {
            _timefile = "timedependence.csv";
        }
        if (options->exists("options.kmcmultiple.carriertype")) {
	    _carriertype = options->get("options.kmcmultiple.carriertype").as<string>();
	}
        else {
	    cout << "WARNING in kmcmultiple: You did not specify a charge carrier type. It will be set to electrons." << endl;
            _carriertype = "e";
        }
        if(_carriertype == "electron" || _carriertype == "electrons" || _carriertype == "e"){
            _carriertype = "e";
            cout << "carrier type: electrons" << endl;
            
        }
        else if(_carriertype == "hole" || _carriertype == "holes" || _carriertype == "h"){
            _carriertype = "h";
            cout << "carrier type: holes" << endl;
            
         }
        else if(_carriertype == "singlet" || _carriertype == "singlets" || _carriertype == "s"){
            _carriertype = "s";
            cout << "carrier type: singlets" << endl;
            
        }
        else if(_carriertype == "triplet" || _carriertype == "triplets" || _carriertype == "t"){
            _carriertype = "t";
            cout << "carrier type: triplets" << endl;
        }
        else {
            _carriertype = "e";
            cout << "Carrier type specification invalid. Setting type to electrons." << endl;
        }
        if (options->exists("options.kmcmultiple.explicitcoulomb")) {
	    _explicitcoulomb = options->get("options.kmcmultiple.explicitcoulomb").as<int>();
	}
        else {
	    cout << "WARNING in kmcmultiple: You did not specify if you want explicit Coulomb interaction to be switched on. It will be switched off." << endl;
            _explicitcoulomb = 0;
        }
        if (options->exists("options.kmcmultiple.temperature")) {
	    _temperature = options->get("options.kmcmultiple.temperature").as<double>();
	}
        else {
	    cout << "WARNING in kmcmultiple: You did not specify a temperature. If no explicit Coulomb interaction is used, this is not a problem, as the rates are read from the state file and temperature is not needed explicitly in the KMC simulation. Otherwise a default value of 300 K is used." << endl;
            _temperature = 300;
        }
        if (options->exists("options.kmcmultiple.rates")) {
	    _rates = options->get("options.kmcmultiple.rates").as<string>();
	}
        else {
	    cout << "Using rates from statefile." << endl;
            _rates = "statefile";
        }
        if(_rates != "statefile" && _rates != "calculate"){
	    cout << "WARNING in kmcmultiple: Invalid option rates. Valid options are 'statefile' or 'calculate'. Setting it to 'statefile'." << endl;
            _rates = "statefile";
        }


        if (options->exists("options.kmcmultiple.tofdirection")) {
	    _tofdirection = options->get("options.kmcmultiple.tofdirection").as<string>();
	}
        if (_tofdirection != "x" && _tofdirection != "y" && _tofdirection != "z"  ) {
	    _tofdirection = "y";
	}
        if (options->exists("options.kmcmultiple.tof")) {
	    _tof = options->get("options.kmcmultiple.tof").as<int>();
	}
        if(_tof == 1){
            cout << "Time of flight experiment: ON" << endl;
            cout << "direction for TOF condition: " << _tofdirection << endl;
        }
        else{
            cout << "Time of flight experiment: OFF (bulk mode)" << endl;
        }
        if (options->exists("options.kmcmultiple.toflength")) {
	    _toflength = options->get("options.kmcmultiple.toflength").as<double>() * 1E-9;
            cout << "Sample length for TOF experiment: " << _toflength*1E9 << " nm" << endl;
	}
        else{
            if(_tof==1){
                cout << "WARNING: time-of-flight mode one, but no sample length (option toflength) defined. Setting it to 50 nm" << endl;
                _toflength = 50.0E-9;
            }
        }
  

        
        

        _filename = filename;
        _outputfile = outputfile;

}


void KMCMultiple::progressbar(double fraction)
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

vector<GNode*> KMCMultiple::LoadGraph()
{
    vector<GNode*> node;
    
    // Load nodes
    votca::tools::Database db;
    db.Open( _filename );
    if(votca::tools::globals::verbose) {cout << "LOADING GRAPH" << endl << "database file: " << _filename << endl; }
    votca::tools::Statement *stmt=NULL;
    if (_carriertype=="h" || _carriertype=="e"){
        stmt= db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnCnN"+_carriertype+", UcNcC"+_carriertype+",eAnion,eNeutral,eCation,UcCnN"+_carriertype+" FROM segments;");
    }
    else if (_carriertype=="s" || _carriertype=="t"){
        stmt= db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnXnN"+_carriertype+", UxNxX"+_carriertype+",eSinglet,eNeutral,eTriplet,UxXnN"+_carriertype+" FROM segments;");
    }
        
    int i=0;
    while (stmt->Step() != SQLITE_DONE)
    {
        GNode *newNode = new GNode();
        node.push_back(newNode);

        int newid = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        node[i]->id = newid;
        vec nodeposition = vec(stmt->Column<double>(2)*1E-9, stmt->Column<double>(3)*1E-9, stmt->Column<double>(4)*1E-9); // converted from nm to m
        node[i]->position = nodeposition;
        node[i]->reorg_intorig = stmt->Column<double>(5); // UnCnN or UnXnN
        node[i]->reorg_intdest = stmt->Column<double>(6); // UcNcC or UxNxX
        double eAnion = stmt->Column<double>(7);
        //double eNeutral = stmt->Column<double>(8);
        double eCation = stmt->Column<double>(9);
        double internalenergy = stmt->Column<double>(10); // UcCnN or UxXnN
        double siteenergy = 0;
        if(_carriertype == "e" || _carriertype=="s")
        {
            siteenergy = eAnion + internalenergy;
        }
        else if(_carriertype == "h" || _carriertype=="t")
        {
            siteenergy = eCation + internalenergy;
        }
        
        node[i]->siteenergy = siteenergy;
        if (votca::tools::wildcmp(_injection_name.c_str(), name.c_str()))
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
    if(votca::tools::globals::verbose) { cout << "segments: " << node.size() << endl; }
    
    // Load pairs and rates
    int numberofpairs = 0;
    stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', rate12"+_carriertype+" AS 'rate', drX, drY, drZ, Jeff2"+_carriertype+", lO"+_carriertype+" FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', rate21"+_carriertype+" AS 'rate', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ', Jeff2"+_carriertype+", lO"+_carriertype+" FROM pairs ORDER BY segment1;");
    while (stmt->Step() != SQLITE_DONE)
    {
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);

        double rate12 = stmt->Column<double>(2);

        vec dr = vec(stmt->Column<double>(3)*1E-9, stmt->Column<double>(4)*1E-9, stmt->Column<double>(5)*1E-9); // converted from nm to m
        double Jeff2 = stmt->Column<double>(6);
        double reorg_out = stmt->Column<double>(7); 
        node[seg1]->AddEvent(seg2,rate12,dr,Jeff2,reorg_out);
        numberofpairs ++;
    }    
    delete stmt;

    if(votca::tools::globals::verbose) { cout << "pairs: " << numberofpairs/2 << endl; }
    
    // Calculate initial escape rates !!!THIS SHOULD BE MOVED SO THAT IT'S NOT DONE TWICE IN CASE OF COULOMB INTERACTION!!!
    for(unsigned int i=0; i<node.size(); i++)
    {
        node[i]->InitEscapeRate();
    }
    return node;
}

KMCMultiple::CoulombMap KMCMultiple::LoadCoulomb(int numberofnodes)
{
    // Load minimum coulomb energy
    votca::tools::Database db;
    db.Open( _filename );
    votca::tools::Statement *stmt1 = db.Prepare("SELECT MIN(coulomb_"+_carriertype+_carriertype+") FROM coulomb;");
    stmt1->Step();
    double mincoulomb = stmt1->Column<double>(0);
    
    cout << "Minimum Coulomb energy of " << mincoulomb << " eV chosen as offset." << endl;
    
    
    // Load Coulomb energies
    CoulombMap coulomb;
    cout << "Loading Coulomb energies from database file " << _filename << endl;
    votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', coulomb_"+_carriertype+_carriertype+" FROM coulomb UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', coulomb_"+_carriertype+_carriertype+" FROM coulomb ORDER BY segment1;");

    int numberofinteractions = 0;
    while (stmt->Step() != SQLITE_DONE)
    {
        
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);
        double coulombenergy = stmt->Column<double>(2) - mincoulomb;
        coulomb[seg1+numberofnodes*seg2] = coulombenergy;
        numberofinteractions ++;
    }
    cout << "    " << numberofinteractions << " Coulomb interactions energies loaded." << endl;
    return coulomb;
}

void ResetForbidden(vector<int> &forbiddenid)
{
    forbiddenid.clear();
}

void AddForbidden(int id, vector<int> &forbiddenid)
{
    forbiddenid.push_back(id);
}

int Forbidden(int id, vector<int> forbiddenlist)
{
    // cout << "forbidden list has " << forbiddenlist.size() << " entries" << endl;
    int forbidden = 0;
    for (unsigned int i=0; i< forbiddenlist.size(); i++)
    {
        if(id == forbiddenlist[i])
        {
            forbidden = 1;
            //cout << "ID " << id << " has been found as element " << i << " (" << forbiddenlist[i]<< ") in the forbidden list." << endl;
            break;
        }
    }
    return forbidden;
}

int Surrounded(GNode* node, vector<int> forbiddendests)
{
    int surrounded = 1;
    for(unsigned int i=0; i<node->event.size(); i++)
    {
        int thisevent_possible = 1;
        for(unsigned int j=0; j<forbiddendests.size(); j++)
        {
            if(node->event[i].destination == forbiddendests[j])
            {
                thisevent_possible = 0;
                break;
            }
        }
        if(thisevent_possible == 1)
        {
            surrounded = 0;
            break;
        }
    }
    return surrounded;
}

void printtime(int seconds_t)
{
    int seconds = seconds_t;
    int minutes = 0;
    int hours = 0;
    while(seconds / 60 >= 1)
    {
        seconds -= 60;
        minutes +=  1;
    }
    while(minutes / 60 >= 1)
    {
        minutes -= 60;
        hours +=  1;
    }
    char buffer [50];
    sprintf(buffer, "%d:%02d:%02d",hours,minutes,seconds);
    printf("%s",buffer);
}


void KMCMultiple::InitBoxSize(vector<GNode*> node)
{
    cout << endl << "Analysing size of the simulation cell." << endl;
    double minX = 0;
    double minY = 0;
    double minZ = 0;
    double maxX = 0;
    double maxY = 0;
    double maxZ = 0;
    double mindX = 777777777777;
    double mindY = 777777777777;
    double mindZ = 777777777777;
    for(unsigned int i = 0; i<node.size(); i++)
    {
        double posX = node[i]->position.x();
        double posY = node[i]->position.y();
        double posZ = node[i]->position.z();
        if(posX < minX) minX = posX;
        if(posY < minY) minY = posY;
        if(posZ < minZ) minZ = posZ;
        if(posX > maxX) maxX = posX;
        if(posY > maxY) maxY = posY;
        if(posZ > maxZ) maxZ = posZ;
        for(unsigned int j = 0; j<node[i]->event.size(); j++)
        {
            if(node[i]->event[j].dr.x() < mindX && node[i]->event[j].dr.x() > 0) mindX = node[i]->event[j].dr.x();
            if(node[i]->event[j].dr.y() < mindY && node[i]->event[j].dr.y() > 0) mindY = node[i]->event[j].dr.y();
            if(node[i]->event[j].dr.z() < mindZ && node[i]->event[j].dr.z() > 0) mindZ = node[i]->event[j].dr.z();
        }   
    }
    _boxsizeX = maxX-minX + mindX;
    _boxsizeY = maxY-minY + mindY;
    _boxsizeZ = maxZ-minZ + mindZ;
    
    double maxdist = _boxsizeX;
    if(_boxsizeY > maxdist) maxdist = _boxsizeY;
    if(_boxsizeZ > maxdist) maxdist = _boxsizeZ;
    maxdist = maxdist /2;
    
    // cout << "lattice constants (X,Y,Z): " << mindX << ", " << mindY << ", " << mindZ << endl;
    cout << "cell dimensions: " << _boxsizeX*1E9 << " x " << _boxsizeY*1E9 << " x " << _boxsizeZ*1E9 << "nm^3" << endl;
    cout << "spatial electron density: " << _numberofcharges/_boxsizeX/_boxsizeY/_boxsizeZ << " m^-3" << endl;
    _coulomboffset = 1.4399644850445791e-09 / maxdist; 
    cout << "Coulomb energy offset: " << _coulomboffset << " eV (for the maximal distance " << maxdist*1E9 << " nm)." << endl;
    
}

void KMCMultiple::InitialRates(vector<GNode*> node)
{
    cout << endl <<"Calculating initial Marcus rates." << endl;
    cout << "    Temperature T = " << _temperature << " K." << endl;
    cout << "    Field: (" << _fieldX << ", " << _fieldY << ", " << _fieldZ << ") V/m" << endl;
    if (_carriertype == "e")
    {
        _q = -1.0;
    }
    else if (_carriertype == "h")
    {
        _q = 1.0;
    }
    else if (_carriertype == "s")
    {
        _q = 0.0;
    }
    else if (_carriertype == "t")
    {
        _q = 0.0;
    }
    cout << "    carriertype: " << _carriertype << " (charge: " << _q << " eV)" << endl;
    int numberofsites = node.size();
    cout << "    Rates for "<< numberofsites << " sites are computed." << endl;
    double maxreldiff = 0;
    int totalnumberofrates = 0;
    for(int i=0; i<numberofsites; i++)
    {
        int numberofneighbours = node[i]->event.size();
        for(int j=0; j<numberofneighbours; j++)
        {
            double dX =  node[i]->event[j].dr.x();
            double dY =  node[i]->event[j].dr.y();
            double dZ =  node[i]->event[j].dr.z();
            double dG_Field = _q * (dX*_fieldX +  dY*_fieldY + dZ*_fieldZ);
            
            double destindex =  node[i]->event[j].destination;
            double reorg = node[i]->reorg_intorig + node[destindex]->reorg_intdest + node[i]->event[j].reorg_out;
            
            double dG_Site = node[destindex]->siteenergy - node[i]->siteenergy;

            double J2 = node[i]->event[j].Jeff2;

            double dG = dG_Site - dG_Field;
            
            double rate = 2*Pi/hbar * J2/sqrt(4*Pi*reorg*kB*_temperature) * exp(-(dG+reorg)*(dG+reorg) / (4*reorg*kB*_temperature));
            
            // calculate relative difference compared to values in the table
            double reldiff = (node[i]->event[j].rate - rate) / node[i]->event[j].rate;
            if (reldiff > maxreldiff) {maxreldiff = reldiff;}
            reldiff = (node[i]->event[j].rate - rate) / rate;
            if (reldiff > maxreldiff) {maxreldiff = reldiff;}
            
            // set rates to calculated values
            node[i]->event[j].rate = rate;
            node[i]->event[j].initialrate = rate;
            
            totalnumberofrates ++;

        }
        
        // Initialise escape rates
        for(unsigned int i=0; i<node.size(); i++)
        {
            node[i]->InitEscapeRate();
        }

    }
    cout << "    " << totalnumberofrates << " rates have been calculated." << endl;
    if(maxreldiff < 0.01)
    {
        cout << "    Good agreement with rates in the state file. Maximal relative difference: " << maxreldiff*100 << " %"<< endl;
    }
    else
    {
        cout << "    WARNING: Rates differ from those in the state file up to " << maxreldiff*100 << " %." << " If the rates in the state file are calculated for a different temperature/field or if they are not Marcus rates, this is fine. Otherwise something might be wrong here." << endl;
    }
}

void KMCMultiple::RateUpdateCoulomb(vector<GNode*> &node,  vector< Chargecarrier* > &carrier, CoulombMap &coulomb)
{
    // Calculate new rates for all possible events, i.e. for the possible hoppings of all occupied nodes
    //cout << "Updating Coulomb interaction part of rates." << endl;
    //#pragma omp parallel for
    for(unsigned int cindex=0; cindex<carrier.size(); cindex++)
    {
        //cout << "  carrier no. "<< cindex+1 << endl;
        GNode *node_i = carrier[cindex]->node;
        double escaperate = 0;
        //#pragma omp parallel for
        for(unsigned int destindex=0; destindex<node_i->event.size(); destindex++)
        {
            int destid = node_i->event[destindex].destination;
            GNode *node_j = node[destid];
            if(node_j->occupied == 1)
            {  // in principal shouldn't be needed:
               node_i->event[destindex].rate = 0;
               // escaperate += node_i->event[destindex].rate;
            }
            else
            {
                //cout << "    event i->j: " << node_i->id+1 << " -> " << node_j->id+1 << endl;
                // getting the correction factor for the rate i->j (node_i -> node_j)
                // now calculating the contribution for all neighbouring charges
                double coulombsum = 0;
                
                int dimension = node.size();
                //#pragma omp parallel for reduction(+:coulombsum)
                for(unsigned int ncindex=0; ncindex<carrier.size(); ncindex++)
                {
                    if(_explicitcoulomb == 2) // raw Coulomb
                    {
                        if(ncindex != cindex)
                        {
                            // - E_ik
                            vec distancevec = carrier[ncindex]->node->position - node_i->position;
                            double epsR = 1;
                            double dX = std::abs(distancevec.x());
                            double dY = std::abs(distancevec.y());
                            double dZ = std::abs(distancevec.z());
                            if (dX > _boxsizeX/2)
                            {
                                dX = _boxsizeX - dX;
                            }
                            if (dY > _boxsizeY/2)
                            {
                                dY = _boxsizeY - dY;
                            }
                            if (dZ > _boxsizeZ/2)
                            {
                                dZ = _boxsizeZ - dZ;
                            }
                            double distance = sqrt(dX*dX+dY*dY+dZ*dZ);
                            double rawcoulombcontribution = 1.4399644850445791e-09 / epsR / distance - _coulomboffset;
                            coulombsum -= rawcoulombcontribution;
                        
                            // + E_jk
                            distancevec = carrier[ncindex]->node->position - node_j->position;
                            dX = std::abs(distancevec.x());
                            dY = std::abs(distancevec.y());
                            dZ = std::abs(distancevec.z());
                            if (dX > _boxsizeX/2)
                            {
                                dX = _boxsizeX - dX;
                            }
                            if (dY > _boxsizeY/2)
                            {
                                dY = _boxsizeY - dY;
                            }
                            if (dZ > _boxsizeZ/2)
                            {
                                dZ = _boxsizeZ - dZ;
                            }
                            distance = sqrt(dX*dX+dY*dY+dZ*dZ);
                            rawcoulombcontribution = 1.4399644850445791e-09 / epsR / distance - _coulomboffset;
                            coulombsum += rawcoulombcontribution;
                        }
                    }
                    
                    else // _explicitcoulomb==1 (partial charges)
                    {
                        CoulombIt coul_iterator;
                        //GNode *node_k = carrier[ncindex]->node;
                        if(ncindex != cindex) // charge doesn't have Coulomb interaction with itself
                        {
                            // - E_ik
                            // check if there is an entry for this interaction
                            unsigned long key = node_i->id + dimension * carrier[ncindex]->node->id;
                            
                            
                            coul_iterator = coulomb.find(key);
                            int additionmade = 0;
                            int substractionmade = 0;
                            double contribution = 0;
                            if( coul_iterator != coulomb.end() )
                            {
                                // if there is an entry, add it to the Coulomb sum
                                //cout << "substracted " << coul_iterator->second << endl;
                                contribution -= coul_iterator->second;
                                substractionmade = 1;
                            }
                            
                            // + E_jk
                            key = node_j->id + dimension * carrier[ncindex]->node->id;
                            coul_iterator = coulomb.find(key);
                            if( coul_iterator != coulomb.end() && substractionmade == 1)
                            {
                                // if there is an entry, add it to the Coulomb sum
                                //cout << "added " << coul_iterator->second << endl;
                                contribution += coul_iterator->second;
                                // cout << " - "<< coul_iterator->second;
                                additionmade = 1;
                            }
                            if(additionmade == 1 && substractionmade == 1)
                            {  // makes sure not to create cutoff-caused unbalanced additions/substractions
                                coulombsum += contribution;
                                // cout << " = "<< contribution << endl;
                            }
                        }
                    // end else mode=partialcharges
                    }

                }
                //double dX =  node_i->event[destindex].dr.x();
                //double dY =  node_i->event[destindex].dr.y();
                //double dZ =  node_i->event[destindex].dr.z();
                //double dG_Field = _q * (dX*_fieldX +  dY*_fieldY + dZ*_fieldZ);
                double reorg = node_i->reorg_intorig + node_j->reorg_intdest + node_i->event[destindex].reorg_out;
                double dG_Site = node_j->siteenergy - node_i->siteenergy;
                //double dG = dG_Site - dG_Field;
                double coulombfactor = exp(-(2*(dG_Site+reorg) * coulombsum + coulombsum*coulombsum) / (4*reorg*kB*_temperature) );
                //if (coulombsum != 0) {
                //cout << "coulombsum = " << coulombsum << endl;
                //cout << "coulombfactor = " << coulombfactor << endl;
                //}
            
                // multiply rates by coulomb factor
                double newrate = node_i->event[destindex].initialrate * coulombfactor;
                // cout << "initial rate: " << node_i->event[destindex].initialrate << endl;
                // cout << "reorg: " << reorg << endl;
                // cout << "dG_Site: " << dG_Site << endl;
                // cout << "dG_Field: " << dG_Field << endl;
                // cout << "coulombsum: " << coulombsum << endl;
                // cout << "coulombfactor: " << coulombfactor << endl;
                // cout << "distance: " << sqrt(dX*dX+dY*dY+dZ*dZ) << endl;
                // cout << "full Coulomb rate: " << newrate << endl << endl;
                // cout << "changed from " << node_i->event[destindex].rate << " to " << newrate << endl;
                node_i->event[destindex].rate = newrate;
            
                escaperate += newrate;
            }
        }
        node_i->escaperate = escaperate;
        
        
    }
}

vector<double> KMCMultiple::RunVSSM(vector<GNode*> node, double runtime, unsigned int numberofcharges, votca::tools::Random2 *RandomVariable, CoulombMap coulomb)
{

    int realtime_start = time(NULL);
    cout << endl << "Algorithm: VSSM for Multiple Charges" << endl;
    cout << "number of charges: " << numberofcharges << endl;
    cout << "number of nodes: " << node.size() << endl;
    string stopcondition;
    unsigned long maxsteps;
    int diffusionsteps = 0;
    unsigned diffusion_stepsize = 10000;
    matrix avgdiffusiontensor;
    avgdiffusiontensor.ZeroMatrix();
    double accumulatedenergy = 0;
    int tdpendencesteps = 1;
    double currentmobility_laststep = 7E77;
    int currentkmcstep_laststep = 0;
    double relaxationtime = 0;
    double relaxationlength = 0;
    int mobilitycheck = 0;
    if(runtime > 100)
    {
        stopcondition = "steps";
        maxsteps = runtime;
        cout << "stop condition: " << maxsteps << " steps." << endl;
    }
    else
    {
        stopcondition = "runtime";
        cout << "stop condition: " << runtime << " seconds runtime." << endl;
        cout << "(If you specify runtimes larger than 100 kmcmultiple assumes that you are specifying the number of steps.)" << endl;
    }
    
    if(numberofcharges > node.size())
    {
        throw runtime_error("ERROR in kmcmultiple: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
    }
    
    fstream traj;
    char trajfile[100];
    strcpy(trajfile, _trajectoryfile.c_str());
    cout << "Writing trajectory to " << trajfile << "." << endl; 
    traj.open (trajfile, fstream::out);
    fstream tfile;
    char timefile[100];
    strcpy(timefile, _timefile.c_str());
    cout << "Writing time dependence of energy and mobility to " << timefile << "." << endl; 
    tfile.open (timefile, fstream::out);
    if(_outputtime != 0)
    {   
        traj << "'time[s]'\t";
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            traj << "'carrier" << i+1 << "_x'\t";    
            traj << "'carrier" << i+1 << "_y'\t";    
            traj << "'carrier" << i+1 << "_z";    
            if(i<numberofcharges-1)
            {
                traj << "'\t";
            }
        }
        traj << endl;
        
        tfile << "'time[s]'\t'energy_per_carrier[eV]\t'mobility[m**2/Vs]'\t'distance_fielddirection[m]'\t'distance_absolute[m]'\t'diffusion_coefficient[m**2]'" << endl;
        
    }
    double outputfrequency = runtime/100;
    double outputstepfrequency = maxsteps/100;
    vector<double> occP(node.size(),0.);
    
    double absolute_field = sqrt(_fieldX*_fieldX + _fieldY*_fieldY + _fieldZ*_fieldZ);

    // Injection
    cout << endl << "injection method: " << _injectionmethod << endl;
    double deltaE = 0;
    double energypercarrier;
    double totalenergy;
    if(_injectionmethod == "equilibrated")
    {
        vector< double > energy;
        for(unsigned int i=0; i<node.size(); i++)
        {
            energy.push_back(node[i]->siteenergy);
        }
        std::sort (energy.begin(), energy.end());
        double fermienergy = 0.5*(energy[_numberofcharges-1]+energy[_numberofcharges]);
                cout << "Energy range in morphology data: [" << energy[0] << ", " << energy[energy.size()-1] << "] eV." << endl;
        cout << "first guess Fermi energy: " << fermienergy << " eV."<< endl;

        cout << "improving guess for fermi energy ";
        
        double probsum;
        double fstepwidth = std::abs(0.01 * fermienergy);
        probsum = 0;
        while(std::abs(probsum - numberofcharges) > 0.1)
        {
            cout << ".";
            totalenergy = 0;
            probsum = 0;
            for(unsigned int i=0; i<energy.size(); i++)
            {
                totalenergy += energy[i] * 1/(exp((energy[i]-fermienergy)/kB/_temperature)+1);
                probsum += 1/(exp((energy[i]-fermienergy)/kB/_temperature)+1);
                //cout << "energy " << energy[i] << " has probability " << 1/(exp((energy[i]-fermienergy)/kB/_temperature)+1) << endl;
            }
            if(std::abs(probsum) > numberofcharges)
            {   // Fermi energy estimate too high
                fermienergy -= fstepwidth;
            }
            else
            {   // Fermi energy estimate too low
                fermienergy += fstepwidth;
            }
            fstepwidth = 0.95 * fstepwidth; // increase precision
        }
        cout << endl << "probsum=" << probsum << "(should be equal to number of charges)." << endl;
        cout << "fermienergy=" << fermienergy << endl;
        cout << "totalenergy=" << totalenergy << endl;
        energypercarrier = totalenergy / probsum; // in theory: probsum=_numberofcharges, but in small systems it can differ significantly
        cout << "Average energy per charge carrier: " << energypercarrier << " eV."<< endl;
        
        double stepwidth = std::abs(fermienergy)/1000;
        int inside = 0;
        while(inside < _numberofcharges)
        {
            inside = 0;
            deltaE += stepwidth;
            for(unsigned int i=0; i<energy.size(); i++)
            {
                if(energy[i] >= energypercarrier-deltaE && energy[i] <= energypercarrier+deltaE)
                {
                    inside += 1;
                }
            }
        }
        cout << inside << " charges in acceptance interval " << energypercarrier << " +/- " << deltaE << "." << endl;
        
    }
    vector< Chargecarrier* > carrier;
    vector<vec> startposition(numberofcharges,vec(0,0,0));
    cout << "looking for injectable nodes..." << endl;
    for (unsigned int i = 0; i < numberofcharges; i++)
    {
        Chargecarrier *newCharge = new Chargecarrier;      
        newCharge->id = i;
        newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];
        int ininterval = 1;
        if (_injectionmethod == "equilibrated") {ininterval = 0;}
        while(newCharge->node->occupied == 1 || newCharge->node->injectable != 1 || ininterval != 1)
        {   // maybe already occupied? or maybe not injectable?
            newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];
            if(_injectionmethod == "equilibrated")
            { // check if charge is in the acceptance interval
                if(newCharge->node->siteenergy >= energypercarrier-0*deltaE && newCharge->node->siteenergy <= energypercarrier+2*deltaE && newCharge->node->occupied == 0 && newCharge->node->injectable == 1)
                {
                    ininterval = 1;
                }
                
            }
        }
        // cout << "selected segment " << newCharge->node->id+1 << " which has energy " << newCharge->node->siteenergy << " within the interval [" << energypercarrier-0*deltaE << ", " << energypercarrier+2*deltaE << "]" << endl;
        newCharge->node->occupied = 1;
        newCharge->dr_travelled = vec(0,0,0);
        newCharge->tof_travelled = 0;
        startposition[i] = newCharge->node->position;
        cout << "starting position for charge " << i+1 << ": segment " << newCharge->node->id+1 << endl;
        carrier.push_back(newCharge);
    }
    
    if(_injectionmethod == "equilibrated")
    {
        cout << "repositioning charges to obtain an equilibrium configuration..." << endl;
        double realisedenergy = 0;
        //while(std::abs(realisedenergy - totalenergy) > 0.01*(std::abs(realisedenergy)+std::abs(totalenergy)))
        //{
            realisedenergy = 0;
            for (unsigned int j = 0; j < numberofcharges; j++)
            {
                realisedenergy += carrier[j]->node->siteenergy;
            }

            for(unsigned int i=0; i<numberofcharges; i++)
            {
                for(unsigned int k=0; k<node.size(); k++)
                {
                    if(std::abs(realisedenergy-carrier[i]->node->siteenergy+node[k]->siteenergy - totalenergy) < std::abs(realisedenergy - totalenergy) && node[k]->occupied == 0)
                    {    //move charge
                         carrier[i]->node->occupied = 0;
                         realisedenergy = realisedenergy-carrier[i]->node->siteenergy+node[k]->siteenergy;
                         carrier[i]->node = node[k];
                         node[k]->occupied = 1;
                    }        
                }
                
            }
        //}    
        cout << "realised energy: " << realisedenergy/numberofcharges << " eV (aimed for " << energypercarrier << " eV)"<<  endl;
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            startposition[i] = carrier[i]->node->position;
            cout << "starting position for charge " << i+1 << ": segment " << carrier[i]->node->id+1 << endl;
        }
    }
    

    double simtime = 0.;
    unsigned long step = 0;
    double nextoutput = outputfrequency;
    unsigned long nextstepoutput = outputstepfrequency;
    double nexttrajoutput = _outputtime;
    unsigned nextdiffstep = diffusion_stepsize;
    
    progressbar(0.);
    vector<int> forbiddennodes;
    vector<int> forbiddendests;
    
    bool timeout = false;
    while(((stopcondition == "runtime" && simtime < runtime) || (stopcondition == "steps" && step < maxsteps)) && timeout == false)
    {
        if((time(NULL) - realtime_start) > _maxrealtime*60.*60.)
        {
            cout  << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime*60*60 +0.5) <<" seconds) has been reached. Stopping here." << endl << endl;
            break;
        }
        double cumulated_rate = 0;
        if(_explicitcoulomb >= 1)
        {
            KMCMultiple::RateUpdateCoulomb(node,carrier,coulomb);
            // cout << "rate11: " << node[0]->event[0].rate << endl;
        }
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            cumulated_rate += carrier[i]->node->getEscapeRate();
        }
        if(cumulated_rate == 0)
        {   // this should not happen: no possible jumps defined for a node
            throw runtime_error("ERROR in kmcmultiple: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
        }
        // go forward in time
        double dt = 0;
        double rand_u = 1-RandomVariable->rand_uniform();
        while(rand_u == 0)
        {
            cout << "WARNING: encountered 0 as a random variable! New try." << endl;
            rand_u = 1-RandomVariable->rand_uniform();
        }
        dt = -1 / cumulated_rate * log(rand_u);
        simtime += dt;
        if(votca::tools::globals::verbose) {cout << "simtime += " << dt << endl << endl;}
        step += 1;
        
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            carrier[i]->node->occupationtime += dt;
        }

        
        ResetForbidden(forbiddennodes);
        int level1step = 0;
        while(level1step == 0 && timeout == false)
        // LEVEL 1
        {
            
            if((time(NULL) - realtime_start) > _maxrealtime*60.*60.)
            {
                cout  << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime*60*60 +0.5) <<" seconds) has been reached while searching escape node (level 1). Stopping here." << endl << endl;
                timeout = true;
                break;
            }

            
            // determine which electron will escape
            GNode* do_oldnode;
            GNode* do_newnode;
            Chargecarrier* do_affectedcarrier;
            
            double u = 1 - RandomVariable->rand_uniform();
            for(unsigned int i=0; i<numberofcharges; i++)
            {
                u -= carrier[i]->node->getEscapeRate()/cumulated_rate;
                if(u <= 0)
                {
                    do_oldnode = carrier[i]->node;
                    do_affectedcarrier = carrier[i];
                    break;
                }
               do_oldnode = carrier[i]->node;
               do_affectedcarrier = carrier[i];
            }
                
            //double maxprob = 0.;
            //double newprob = 0.;
            vec dr;
            if(votca::tools::globals::verbose) {cout << "Charge number " << do_affectedcarrier->id+1 << " which is sitting on segment " << do_oldnode->id+1 << " will escape!" << endl ;}
            if(Forbidden(do_oldnode->id, forbiddennodes) == 1) {continue;}
            
            // determine where it will jump to
            ResetForbidden(forbiddendests);
            while(timeout == false)
            {
            // LEVEL 2
                if(votca::tools::globals::verbose) {cout << "There are " << do_oldnode->event.size() << " possible jumps for this charge:"; }
                
                if((time(NULL) - realtime_start) > _maxrealtime*60.*60.)
                {
                    cout  << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime*60*60 +0.5) <<" seconds) has been reached while searching destination node (level 2). Stopping here." << endl << endl;
                    timeout = true;
                    break;
                }


                do_newnode = NULL;
                u = 1 - RandomVariable->rand_uniform();
                for(unsigned int j=0; j<do_oldnode->event.size(); j++)
                {
                    if(votca::tools::globals::verbose) { cout << " " << do_oldnode->event[j].destination+1 ; }
                    u -= do_oldnode->event[j].rate/do_oldnode->getEscapeRate();
                    if(u <= 0)
                    {
                        do_newnode = node[do_oldnode->event[j].destination];
                        dr = do_oldnode->event[j].dr;
                        break;
                    }
                    do_newnode = node[do_oldnode->event[j].destination];
                    dr = do_oldnode->event[j].dr;
                }

                if(do_newnode == NULL)
                {
                    if(votca::tools::globals::verbose) {cout << endl << "Node " << do_oldnode->id+1  << " is SURROUNDED by forbidden destinations and zero rates. Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl; }
                    AddForbidden(do_oldnode->id, forbiddennodes);
                    break; // select new escape node (ends level 2 but without setting level1step to 1)
                }
                if(votca::tools::globals::verbose) {cout << endl << "Selected jump: " << do_newnode->id+1 << endl; }
                
                // check after the event if this was allowed
                if(Forbidden(do_newnode->id, forbiddendests) == 1)
                {
                    if(votca::tools::globals::verbose) {cout << "Node " << do_newnode->id+1  << " is FORBIDDEN. Now selection new hopping destination." << endl; }
                    continue;
                }

                // if the new segment is unoccupied: jump; if not: add to forbidden list and choose new hopping destination
                if(do_newnode->occupied == 1)
                {
                    if(Surrounded(do_oldnode, forbiddendests) == 1)
                    {
                        if(votca::tools::globals::verbose) {cout << "Node " << do_oldnode->id+1  << " is SURROUNDED by forbidden destinations. Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl; }
                        AddForbidden(do_oldnode->id, forbiddennodes);
                        break; // select new escape node (ends level 2 but without setting level1step to 1)
                    }
                    if(votca::tools::globals::verbose) {cout << "Selected segment: " << do_newnode->id+1 << " is already OCCUPIED. Added to forbidden list." << endl << endl;}
                    AddForbidden(do_newnode->id, forbiddendests);
                    if(votca::tools::globals::verbose) {cout << "Now choosing different hopping destination." << endl; }
                    continue; // select new destination
                }
                else
                {
                    do_newnode->occupied = 1;
                    do_oldnode->occupied = 0;
                    do_affectedcarrier->node = do_newnode;
                    do_affectedcarrier->dr_travelled += dr;
                    level1step = 1;
                    if(votca::tools::globals::verbose) {cout << "Charge has jumped to segment: " << do_newnode->id+1 << "." << endl;}
                    
                    if (_tof == 1)
                    {   // time-of-flight experiment
                        
                        // keep track of travel in TOF direction
                        double dr_tof = 0;
                        if(_tofdirection == "x"){
                            dr_tof = dr.getX();
                        }
                        else if(_tofdirection == "y"){
                            dr_tof = dr.getY();
                        }
                        else {
                            dr_tof = dr.getZ();
                        }
                        
                        do_affectedcarrier->tof_travelled += dr_tof;
                        
                        // boundary crossing
                        if(do_affectedcarrier->tof_travelled > _toflength){
                            cout << "\nTOF length has been crossed in "+_tofdirection+" direction." << endl;
                            cout << "LAST POSITION: " << (do_affectedcarrier->dr_travelled-dr) *1E9 << ", NEW POSITION: " << do_affectedcarrier->dr_travelled *1E9 << endl ; 

                            // now re-inject the carrier to a new random position
                            GNode *temp_node = new GNode;
                            temp_node = do_affectedcarrier->node;
                            do_affectedcarrier->node = node[RandomVariable->rand_uniform_int(node.size())];
                            while(do_affectedcarrier->node->occupied == 1 || do_affectedcarrier->node->injectable != 1)
                            {   // maybe already occupied? or maybe not injectable?
                                do_affectedcarrier->node = node[RandomVariable->rand_uniform_int(node.size())];
                            }
                            do_affectedcarrier->node->occupied = 1;
                            temp_node->occupied = 0;
                            cout << "carrier " << do_affectedcarrier->id << " has been instanteneously moved from node " << temp_node->id+1 << " to node " << do_affectedcarrier->node->id+1 << endl << endl; 
                            
                            // reset TOF length counter
                            do_affectedcarrier->tof_travelled -= _toflength;
                        }
                        
                    }
                    
                    break; // this ends LEVEL 2 , so that the time is updated and the next MC step started
                }

                if(votca::tools::globals::verbose) {cout << "." << endl;}
            // END LEVEL 2
            }
        // END LEVEL 1
        }    
        
        
        if(step > nextdiffstep)       
        {
            nextdiffstep += diffusion_stepsize;
            for(unsigned int i=0; i<numberofcharges; i++) 
            {
                diffusionsteps  += 1;
                avgdiffusiontensor += (startposition[i]+carrier[i]->dr_travelled)|(startposition[i]+carrier[i]->dr_travelled);
            }
        }
        
        if(_outputtime != 0 && simtime > nexttrajoutput)       
        {
            // write to trajectory file
            nexttrajoutput = simtime + _outputtime;
            traj << simtime << "\t";
            for(unsigned int i=0; i<numberofcharges; i++) 
            {
                traj << startposition[i].getX() + carrier[i]->dr_travelled.getX() << "\t";
                traj << startposition[i].getY() + carrier[i]->dr_travelled.getY() << "\t";
                traj << startposition[i].getZ() + carrier[i]->dr_travelled.getZ();
                if (i<numberofcharges-1) 
                {
                    traj << "\t";
                }
                else
                {
                    traj << endl;
                }
            }
            
            // write to time development file
            // a) energy
            double currentenergy = 0;
            // b) mobility
            double currentmobility = 0;
            vec dr_travelled_current = vec (0,0,0);
            if(absolute_field != 0)
            {
                vec avgvelocity_current = vec(0,0,0);
                for(unsigned int i=0; i<numberofcharges; i++)
                {
                    dr_travelled_current += carrier[i]->dr_travelled;
                    accumulatedenergy += carrier[i]->node->siteenergy;
                }
                dr_travelled_current /= numberofcharges;
                currentenergy = accumulatedenergy /tdpendencesteps/numberofcharges;
                avgvelocity_current = dr_travelled_current/simtime; 
                currentmobility = (avgvelocity_current*_field) /absolute_field/absolute_field;
            }
            tfile << simtime << "\t" << currentenergy << "\t" << currentmobility << "\t" << (dr_travelled_current*_field)/absolute_field << "\t" << abs(dr_travelled_current) << endl;
            
            if(currentmobility>0)
            {
                if(std::abs((currentmobility - currentmobility_laststep)/currentmobility) < 0.001)
                {
                    mobilitycheck += step - currentkmcstep_laststep;
                }
                else
                {
                    mobilitycheck = 0;
                }
            }
            if(mobilitycheck >= 1E8 && relaxationtime == 0)
            {
                cout << "\nKMC probably converged at t= " << simtime << endl;
                cout << "    (For the last 10^8 KMC steps the relative difference in mobility was smaller than 0.001.)" << endl;
                relaxationtime = simtime;
                for(unsigned int i=0; i<numberofcharges; i++)
                {
                    double thisX = carrier[i]->dr_travelled.getX();
                    double thisY = carrier[i]->dr_travelled.getY();
                    double thisZ = carrier[i]->dr_travelled.getZ();
                    relaxationlength += sqrt(thisX*thisX+thisY*thisY+thisZ*thisZ);
                }
                relaxationlength /= numberofcharges;
                cout << "    relaxation time: " << relaxationtime << " s." << endl;
                cout << "    relaxation length: " << relaxationlength << " m." << endl;
            }
            
            
            currentmobility_laststep = currentmobility;
            tdpendencesteps += 1;
            currentkmcstep_laststep = step;
            
        }
        if(stopcondition == "runtime" && simtime > nextoutput)
        {
            nextoutput = simtime + outputfrequency;
            progressbar(simtime/runtime);
            cout << " remaining: ";
            printtime(int((runtime/simtime-1) * (int(time(NULL)) - realtime_start))); 
        }
        else if(stopcondition == "steps" && step > nextstepoutput)
        {
            nextstepoutput = step + outputstepfrequency;
            progressbar(double(step)/double(maxsteps));
            cout << " remaining: ";
            printtime(int((double(maxsteps)/double(step)-1) * (int(time(NULL)) - realtime_start))); 
        }
    }
    progressbar(1.);
    
    if(_outputtime != 0)
    {   
        traj.close();
        tfile.close();
    }


    // calculate occupation probabilities from occupation times    
    for(unsigned int j=0; j<node.size(); j++)
    {   
        occP[j] = node[j]->occupationtime / simtime;
    }
    

    cout << endl << "finished KMC simulation after " << step << " steps." << endl;
    cout << "simulated time " << simtime << " seconds." << endl;
    cout << "runtime: ";
    printtime(time(NULL) - realtime_start); 
    cout << endl << endl;
    vec dr_travelled = vec (0,0,0);
    vec avgvelocity = vec(0,0,0);
    for(unsigned int i=0; i<numberofcharges; i++)
    {
        //cout << std::scientific << "    charge " << i+1 << ": " << carrier[i]->dr_travelled/simtime*1e-9 << endl;
        cout << std::scientific << "    charge " << i+1 << ": " << carrier[i]->dr_travelled/simtime << endl;
        dr_travelled += carrier[i]->dr_travelled;
    }
    dr_travelled /= numberofcharges;
    //cout << std::scientific << "  Overall average velocity (m/s): " << dr_travelled/simtime*1e-9 << endl;
    avgvelocity = dr_travelled/simtime; 
    cout << std::scientific << "  Overall average velocity (m/s): " << avgvelocity << endl;

    cout << endl << "Distances travelled (m): " << endl;
    for(unsigned int i=0; i<numberofcharges; i++)
    {
        cout << std::scientific << "    charge " << i+1 << ": " << carrier[i]->dr_travelled << endl;
    }
    
    // calculate mobilities
    double average_mobility = 0;
    if (absolute_field != 0)
    {
        cout << endl << "Mobilities (m^2/Vs): " << endl;
        for(unsigned int i=0; i<numberofcharges; i++)
        {
            //vec velocity = carrier[i]->dr_travelled/simtime*1e-9;
            vec velocity = carrier[i]->dr_travelled/simtime;
            //double absolute_velocity = sqrt(velocity.x()*velocity.x() + velocity.y()*velocity.y() + velocity.z()*velocity.z());
            //cout << std::scientific << "    charge " << i+1 << ": mu=" << absolute_velocity/absolute_field*1E4 << endl;
            cout << std::scientific << "    charge " << i+1 << ": mu=" << (velocity*_field)/absolute_field/absolute_field << endl;
            average_mobility += (velocity*_field) /absolute_field/absolute_field;
        }
        average_mobility /= numberofcharges;
        cout << std::scientific << "  Overall average mobility in field direction <mu>=" << average_mobility << " m^2/Vs (= " << average_mobility*1E4 << "cm^2/Vs)" << endl;
      }
    cout << endl;
    
    // calculate diffusion tensor
    avgdiffusiontensor /= (diffusionsteps*2*simtime*numberofcharges);
    cout<<endl<<"Diffusion tensor averaged over all carriers (m^2/s):" << endl << avgdiffusiontensor << endl;

    matrix::eigensystem_t diff_tensor_eigensystem;
    cout<<endl<<"Eigenvalues: "<<endl<<endl;
    avgdiffusiontensor.SolveEigensystem(diff_tensor_eigensystem);
    for(int i=0; i<=2; i++)
    {
        cout<<"Eigenvalue: "<<diff_tensor_eigensystem.eigenvalues[i]<<endl<<"Eigenvector: ";
               
        cout<<diff_tensor_eigensystem.eigenvecs[i].x()<<"   ";
        cout<<diff_tensor_eigensystem.eigenvecs[i].y()<<"   ";
        cout<<diff_tensor_eigensystem.eigenvecs[i].z()<<endl<<endl;
    }
    
    // calculate average mobility from the Einstein relation
    if (absolute_field == 0)
    {
       cout << "The following value is calculated using the Einstein relation and assuming an isotropic medium" << endl;
       double avgD  = 1./3. * (diff_tensor_eigensystem.eigenvalues[0] + diff_tensor_eigensystem.eigenvalues[1] + diff_tensor_eigensystem.eigenvalues[2] );
       average_mobility = std::abs(avgD / kB / _temperature);
       cout << std::scientific << "  Overall average mobility <mu>=" << average_mobility << " m^2/Vs (= " << average_mobility*1E4 << "cm^2/Vs)\n\n" << endl;
    }
    
    
    // check if Einstein relation is fulfilled
    string direction = "";
    double field = 0;
    if(_fieldX != 0 && _fieldY == 0 && _fieldZ == 0) {direction = "x"; field = _fieldX;}
    else if(_fieldX == 0 && _fieldY != 0 && _fieldZ == 0) {direction = "y"; field = _fieldY;}
    else if(_fieldX == 0 && _fieldY == 0 && _fieldZ != 0) {direction = "z"; field = _fieldZ;}
    if(direction != "")
    {
        cout << "components of the mobility tensor in " << direction << " direction (m^2/Vs):" << endl;
        double mu1 = avgvelocity.getX()/field;
        double mu2 = avgvelocity.getY()/field;
        double mu3 = avgvelocity.getZ()/field;
        cout << "mu_x" << direction << " = " << mu1 << endl;
        cout << "mu_y" << direction << " = " << mu2<< endl;
        cout << "mu_z" << direction << " = " << mu3 << endl;
        
        //cout << "\nideality factor for the Einstein relation in " << direction << " direction." << endl;
        //double D1;
        //double D2;
        //double D3;
        //if(direction == "x"){D1 = avgdiffusiontensor.get(0,0);D2 = avgdiffusiontensor.get(1,0);D3 = avgdiffusiontensor.get(2,0);}
        //else if(direction == "y"){D1 = avgdiffusiontensor.get(0,1);D2 = avgdiffusiontensor.get(1,1);D3 = avgdiffusiontensor.get(2,1);}
        //else if(direction == "z"){D1 = avgdiffusiontensor.get(0,2);D2 = avgdiffusiontensor.get(1,2);D3 = avgdiffusiontensor.get(2,2);}
        //cout << "g_x" << direction << " = " << (D1/mu1) / (kB*_temperature) << endl;
        //cout << "g_y" << direction << " = " << (D2/mu2) / (kB*_temperature) << endl;
        //cout << "g_z" << direction << " = " << (D3/mu3) / (kB*_temperature) << endl;
    }
    
    if(_outputtime != 0 && relaxationtime != 0)
    {
        cout << "\nKMC probably converged at t= " << relaxationtime << endl;
        cout << "    (For the last 10^8 KMC steps the relative difference in mobility was smaller than 0.001.)" << endl;
        cout << "    convergence time: " << relaxationtime << " s." << endl;
        cout << "    convergence length: " << relaxationlength << " m." << endl;
    }
    else if(_outputtime != 0 && relaxationtime == 0)
    {
        //cout << "\nKMC probably has not converged yet." << endl;
    }


    
    return occP;
}


void KMCMultiple::WriteOcc(vector<double> occP, vector<GNode*> node)
{
    votca::tools::Database db;
    cout << "Opening for writing " << _filename << endl;
	db.Open(_filename);
	db.Exec("BEGIN;");
	votca::tools::Statement *stmt = db.Prepare("UPDATE segments SET occP"+_carriertype+" = ? WHERE _id = ?;");
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
    std::cout << "-----------------------------------" << std::endl;      
    std::cout << "      KMC FOR MULTIPLE CHARGES" << std::endl;
    std::cout << "-----------------------------------" << std::endl << std::endl;      
 
    // Initialise random number generator
    if(votca::tools::globals::verbose) { cout << endl << "Initialising random number generator" << endl; }
    srand(_seed); // srand expects any integer in order to initialise the random number generator
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());
    
    vector<GNode*> node;
    node = KMCMultiple::LoadGraph();
    CoulombMap coulomb;
    if(_explicitcoulomb == 1)
    {
        cout << endl << "Explicit Coulomb Interaction: ON (Partial Charges)." << endl << "[explicitcoulomb=1]" << endl;
        KMCMultiple::InitialRates(node);
        coulomb = KMCMultiple::LoadCoulomb(node.size());
    }
    else if(_explicitcoulomb == 2)
    {
        cout << endl << "Explicit Coulomb Interaction: ON (Raw Coulomb Interaction)." << endl << "[explicitcoulomb=2]" << endl;
        KMCMultiple::InitialRates(node);
        KMCMultiple::InitBoxSize(node);
    }
    else
    {
        cout << endl << "Explicit Coulomb Interaction: OFF." << endl;
        if(_rates == "calculate")
        {
            cout << "Calculating rates (i.e. rates from state file are not used)." << endl;
            KMCMultiple::InitialRates(node);
        }
        else
        {
            cout << "Using rates from state file." << endl;
        }
    }
    vector<double> occP(node.size(),0.);

    occP = KMCMultiple::RunVSSM(node, _runtime, _numberofcharges, RandomVariable, coulomb);
    
    // write occupation probabilites
    KMCMultiple::WriteOcc(occP, node);
    
    return true;
}

}}


#endif	/* __VOTCA_KMC_MULTIPLE_H */