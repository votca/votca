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

#include "calculators/kmclifetime.h"
#include <votca/tools/property.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>

using namespace std;

namespace votca {
    namespace xtp {

        namespace conv=votca::tools::conv;
        const double Pi = boost::math::constants::pi<double>();
        
        
        void KMCLifetime::Initialize(const char *filename, Property *options, const char *outputfile) {
            if (options->exists("options.kmclifetime.insertions")) {
                _insertions = options->get("options.kmclifetime.insertions").as<unsigned int>();
            } else {
                throw runtime_error("Error in kmclifetime: total number of insertions not provided");
            }
            if (options->exists("options.kmcmultiple.maxrealtime")) {
	    _maxrealtime = options->get("options.kmcmultiple.maxrealtime").as<double>();
            }
            else{
            _maxrealtime = 1E10; // maximal real time in hours
            }
            if (options->exists("options.kmclifetime.seed")) {
                _seed = options->get("options.kmclifetime.seed").as<int>();
            } else {
                throw runtime_error("Error in kmclifetime: seed is not provided");
            }

            if (options->exists("options.kmclifetime.numberofcharges")) {
                _numberofcharges = options->get("options.kmclifetime.numberofcharges").as<int>();
            } else {
                throw runtime_error("Error in kmclifetime: number of charges is not provided");
            }

            if (options->exists("options.kmclifetime.injection")) {
                _injection_name = options->get("options.kmclifetime.injection").as<string>();
            } else {
                throw runtime_error("Error in kmclifetime: injection pattern is not provided");
            }
            if (options->exists("options.kmclifetime.lifetime")) {
                _lifetimefile = options->get("options.kmclifetime.lifetime").as<string>();
            } else {
                throw runtime_error("Error in kmclifetime: Lifetime file not provided");
            }

            if (options->exists("options.kmclifetime.trajectoryfile")) {
                _trajectoryfile = options->get("options.kmclifetime.trajectoryfile").as<string>();
            } else {
                _trajectoryfile = "trajectory.csv";
            }
            if (_trajectoryfile == "") {
                _trajectoryfile = "trajectory.csv";
            }
            if (options->exists("options.kmclifetime.carrierenergy")) {
                _do_carrierenergy=options->get("options.kmclifetime.carrierenergy.run").as<bool>();
                _energy_outputfile = options->get("options.kmclifetime.carrierenergy.outputfile").as<string>();
                _alpha = options->get("options.kmclifetime.carrierenergy.alpha").as<double>();
                _outputsteps = options->get("options.kmclifetime.carrierenergy.outputsteps").as<unsigned>();
                
            } else {
                _do_carrierenergy=false;
            }
            
          
            if (options->exists("options.kmclifetime.carriertype")) {
                _carriertype = options->get("options.kmclifetime.carriertype").as<string>();
            } else {
                cout << "WARNING in kmclifetime: You did not specify a charge carrier type. It will be set to singlets." << endl;
                _carriertype = "s";
            }
            if (_carriertype == "electron" || _carriertype == "electrons" || _carriertype == "e") {
                _carriertype = "e";
                cout << "carrier type: electrons" << endl;

            } else if (_carriertype == "hole" || _carriertype == "holes" || _carriertype == "h") {
                _carriertype = "h";
                cout << "carrier type: holes" << endl;

            } else if (_carriertype == "singlet" || _carriertype == "singlets" || _carriertype == "s") {
                _carriertype = "s";
                cout << "carrier type: singlets" << endl;

            } else if (_carriertype == "triplet" || _carriertype == "triplets" || _carriertype == "t") {
                _carriertype = "t";
                cout << "carrier type: triplets" << endl;
            } else {
                _carriertype = "e";
                cout << "Carrier type specification invalid. Setting type to electrons." << endl;
            }

            if (options->exists("options.kmclifetime.temperature")) {
                _temperature = options->get("options.kmclifetime.temperature").as<double>();
            } else {
                cout << "WARNING in kmclifetime: You did not specify a temperature. A default value of 300 K is used." << endl;     
                _temperature = 300;
            }
            if (options->exists("options.kmclifetime.rates")) {
                _rates = options->get("options.kmclifetime.rates").as<string>();
            } else {
                cout << "Using rates from statefile." << endl;
                _rates = "statefile";
            }
            if (_rates != "statefile" && _rates != "calculate") {
                cout << "WARNING in kmclifetime: Invalid option rates. Valid options are 'statefile' or 'calculate'. Setting it to 'statefile'." << endl;
                _rates = "statefile";
            }


            _filename = filename;
            _outputfile = outputfile;

        }




        vector<GNode*> KMCLifetime::LoadGraph() {
            vector<GNode*> node;

            // Load nodes
            votca::tools::Database db;
            db.Open(_filename);
           
            votca::tools::Statement *stmt;
            if (_carriertype == "h" || _carriertype == "e") {
                stmt = db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnCnN" + _carriertype + ", UcNcC" + _carriertype + ",eAnion,eNeutral,eCation,UcCnN" + _carriertype + " FROM segments;");
            } else if (_carriertype == "s" || _carriertype == "t") {
                stmt = db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnXnN" + _carriertype + ", UxNxX" + _carriertype + ",eSinglet,eNeutral,eTriplet,UxXnN" + _carriertype + " FROM segments;");
            } else{
                throw runtime_error("Carriertype "+_carriertype+" not known.");
            }

            int i = 0;
            while (stmt->Step() != SQLITE_DONE) {
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
                if (_carriertype == "e" || _carriertype == "s") {
                    siteenergy = eAnion + internalenergy;
                } else if (_carriertype == "h" || _carriertype == "t") {
                    siteenergy = eCation + internalenergy;
                }

                node[i]->siteenergy = siteenergy;
                if (votca::tools::wildcmp(_injection_name.c_str(), name.c_str())) {
                    node[i]->injectable = 1;
                } else {
                    node[i]->injectable = 0;
                }
                i++;
            }
            delete stmt;
            
            ReadLifetimeFile(_lifetimefile,node);

            // Load pairs and rates
            int numberofpairs = 0;
            stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', rate12" + _carriertype + " AS 'rate', drX, drY, drZ, Jeff2" + _carriertype + ", lO" + _carriertype + " FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', rate21" + _carriertype + " AS 'rate', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ', Jeff2" + _carriertype + ", lO" + _carriertype + " FROM pairs ORDER BY segment1;");
            while (stmt->Step() != SQLITE_DONE) {
                int seg1 = stmt->Column<int>(0);
                int seg2 = stmt->Column<int>(1);

                double rate12 = stmt->Column<double>(2);

                vec dr = vec(stmt->Column<double>(3)*1E-9, stmt->Column<double>(4)*1E-9, stmt->Column<double>(5)*1E-9); // converted from nm to m
                double Jeff2 = stmt->Column<double>(6);
                double reorg_out = stmt->Column<double>(7);
                node[seg1]->AddEvent(seg2, rate12, dr, Jeff2, reorg_out);
                numberofpairs++;
            }
            delete stmt;

            if (votca::tools::globals::verbose) {
                cout << "pairs: " << numberofpairs / 2 << endl;
            }

            //get boxsize
            double boxsize = 1.0;
            stmt = db.Prepare("SELECT box11,box12,box13,box21,box22,box23,box31,box32,box33 FROM frames");
            while (stmt->Step() != SQLITE_DONE) {
                vec vecx = vec(stmt->Column<double>(0), stmt->Column<double>(1), stmt->Column<double>(2));
                vec vecy = vec(stmt->Column<double>(3), stmt->Column<double>(4), stmt->Column<double>(5));
                vec vecz = vec(stmt->Column<double>(6), stmt->Column<double>(7), stmt->Column<double>(8));
                boxsize = std::abs(((vecx^vecy) * vecz));
            }

            cout << "spatial density: " << _numberofcharges / boxsize << " nm^-3" << endl;

            for (unsigned int i = 0; i < node.size(); i++) {
                node[i]->InitEscapeRate();
            }
            return node;
        }

        void KMCLifetime::ResetForbiddenlist(vector<int> &forbiddenid) {
            forbiddenid.clear();
        }

        void KMCLifetime::AddtoForbiddenlist(int id, vector<int> &forbiddenid) {
            forbiddenid.push_back(id);
        }

        bool KMCLifetime::CheckForbidden(int id,const vector<int> &forbiddenlist) {
            // cout << "forbidden list has " << forbiddenlist.size() << " entries" << endl;
            bool forbidden = false;
            for (unsigned int i = 0; i < forbiddenlist.size(); i++) {
                if (id == forbiddenlist[i]) {
                    forbidden = true;
                    //cout << "ID " << id << " has been found as element " << i << " (" << forbiddenlist[i]<< ") in the forbidden list." << endl;
                    break;
                }
            }
            return forbidden;
        }

        bool KMCLifetime::CheckSurrounded(GNode* node,const vector<int> & forbiddendests) {
            bool surrounded = true;
            for (unsigned  i = 0; i < node->event.size(); i++) {
                bool thisevent_possible = true;
                for (unsigned int j = 0; j < forbiddendests.size(); j++) {
                    if (node->event[i].destination == forbiddendests[j]) {
                        thisevent_possible = false;
                        break;
                    }
                }
                if (thisevent_possible == true) {
                    surrounded = false;
                    break;
                }
            }
            return surrounded;
        }

        void KMCLifetime::printtime(int seconds_t) {
            int seconds = seconds_t;
            int minutes = 0;
            int hours = 0;
            while (seconds / 60 >= 1) {
                seconds -= 60;
                minutes += 1;
            }
            while (minutes / 60 >= 1) {
                minutes -= 60;
                hours += 1;
            }
            char buffer [50];
            sprintf(buffer, "%d:%02d:%02d", hours, minutes, seconds);
            printf("%s", buffer);
        }
        
        void KMCLifetime::ReadLifetimeFile(string filename,vector<GNode*> node){
            Property xml;
            load_property_from_xml(xml, filename);
            list<Property*> jobProps = xml.Select("lifetimes.site");
            if (jobProps.size()!=node.size()){
                throw  runtime_error((boost::format("The number of sites in the sqlfile: %i does not match the number in the lifetimefile: %i") % node.size() % jobProps.size()).str());
                
            }
            
            for (list<Property*> ::iterator  it = jobProps.begin(); it != jobProps.end(); ++it) {
                int site_id =(*it)->getAttribute<int>("id");
                double lifetime=boost::lexical_cast<double>((*it)->value());
                bool check=false;
                for (unsigned i=0;i<node.size();i++){
                    if (node[i]->id==site_id-1 && !(node[i]->hasdecay)){
                        node[i]->AddDecayEvent(1.0/lifetime);
                        check=true;
                        break;
                    }
                    else if(node[i]->id==site_id && node[i]->hasdecay){
                        throw runtime_error((boost::format("Node %i appears twice in your list") %site_id).str());
                    } 
                    
                }
                if (!check){
                throw runtime_error((boost::format("Site from file with id: %i not found in sql") %site_id).str());
                }
            }
            
        }

        void KMCLifetime::InitialRates(vector<GNode*> node) {
            cout << endl << "Calculating initial Marcus rates." << endl;
            cout << "    Temperature T = " << _temperature << " K." << endl;
           
            cout << "    carriertype: " << _carriertype << endl;
            unsigned numberofsites = node.size();
            cout << "    Rates for " << numberofsites << " sites are computed." << endl;
            double maxreldiff = 0;
            int totalnumberofrates = 0;
            for (unsigned int i = 0; i < numberofsites; i++) {
                unsigned numberofneighbours = node[i]->event.size();
                for (unsigned int j = 0; j < numberofneighbours; j++) {
                    if(node[i]->event[j].decayevent){
                        //if event is a decay event there is no point in calculating its rate, because it already has that from the reading in.
                        continue;
                    }
                   
                    

                    double destindex = node[i]->event[j].destination;
                    double reorg = node[i]->reorg_intorig + node[destindex]->reorg_intdest + node[i]->event[j].reorg_out;

                    double dG = node[destindex]->siteenergy - node[i]->siteenergy;

                    double J2 = node[i]->event[j].Jeff2;

                    double rate = 2 * Pi / conv::hbar * J2 / sqrt(4 * Pi * reorg * conv::kB * _temperature) * exp(-(dG + reorg)*(dG + reorg) / (4 * reorg * conv::kB * _temperature));

                    // calculate relative difference compared to values in the table
                    double reldiff = (node[i]->event[j].rate - rate) / node[i]->event[j].rate;
                    if (reldiff > maxreldiff) {
                        maxreldiff = reldiff;
                    }
                    reldiff = (node[i]->event[j].rate - rate) / rate;
                    if (reldiff > maxreldiff) {
                        maxreldiff = reldiff;
                    }

                    // set rates to calculated values
                    node[i]->event[j].rate = rate;
                    node[i]->event[j].initialrate = rate;

                    totalnumberofrates++;

                }

                // Initialise escape rates
                for (unsigned int i = 0; i < node.size(); i++) {
                    node[i]->InitEscapeRate();
                }

            }
            cout << "    " << totalnumberofrates << " rates have been calculated." << endl;
            if (maxreldiff < 0.01) {
                cout << "    Good agreement with rates in the state file. Maximal relative difference: " << maxreldiff * 100 << " %" << endl;
            } else {
                cout << "    WARNING: Rates differ from those in the state file up to " << maxreldiff * 100 << " %." << " If the rates in the state file are calculated for a different temperature/field or if they are not Marcus rates, this is fine. Otherwise something might be wrong here." << endl;
            }
            return;
        }
        
        double KMCLifetime::Promotetime(double cumulated_rate,votca::tools::Random2 * RandomVariable){
            double dt = 0;
                double rand_u = 1 - RandomVariable->rand_uniform();
                while (rand_u == 0) {
                    cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                    rand_u = 1 - RandomVariable->rand_uniform();
                }
                dt = -1 / cumulated_rate * log(rand_u);
            return dt;
        }
        
        GLink* KMCLifetime::ChooseHoppingDest(GNode* node,votca::tools::Random2 * RandomVariable){
            double u = 1 - RandomVariable->rand_uniform();
            
            for (unsigned int j = 0; j < node->event.size(); j++) {
                u -= node->event[j].rate / node->getEscapeRate();
                if (u <= 0 || j==node->event.size()-1) {                   
                    return &(node->event[j]);
                }
            }
            return NULL;
        }
        

        vector<double>  KMCLifetime::RunVSSM(vector<GNode*> node, unsigned int insertions, unsigned int numberofcharges, votca::tools::Random2 * RandomVariable) {

            int realtime_start = time(NULL);
            cout << endl << "Algorithm: VSSM for Multiple Charges with finite Lifetime" << endl;
            cout << "number of charges: " << numberofcharges << endl;
            cout << "number of nodes: " << node.size() << endl;
            
            
            

            if (numberofcharges > node.size()) {
                throw runtime_error("ERROR in kmclifetime: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
            }

            vector<double> occP(node.size(), 0.);

            fstream traj;
            fstream energyfile;

            cout << "Writing trajectory to " <<  _trajectoryfile << "." << endl; 
            traj.open ( _trajectoryfile.c_str(), fstream::out);
            traj << "#Simtime [s]\t Insertion\t Carrier ID\t Lifetime[s]\tSteps\t Last Segment\t x_travelled[m]\t y_travelled[m]\t z_travelled[m]"<<endl;
            
            if(_do_carrierenergy){

                cout << "Tracking the energy of one charge carrier and exponential average with alpha=" << _alpha << " to "<<_energy_outputfile << endl;
                energyfile.open(_energy_outputfile.c_str(),fstream::out);
                energyfile << "Simtime [s]\tSteps\tCarrier ID\tEnergy_a="<<_alpha<<"[nm]"<<endl;
            }
            
            
            
            // Injection
            cout << endl << "injection method: " << _injectionmethod << endl;
            
            vector< Chargecarrier* > carrier;
            vector<vec> startposition(numberofcharges, vec(0, 0, 0));
            cout << "looking for injectable nodes..." << endl;
            for (unsigned int i = 0; i < numberofcharges; i++) {
                Chargecarrier *newCharge = new Chargecarrier;
                newCharge->id = i;
                do {
                newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];     
                }
                while (newCharge->node->occupied == 1 || newCharge->node->injectable != 1 ); // maybe already occupied? or maybe not injectable?
                    
                
                // cout << "selected segment " << newCharge->node->id+1 << " which has energy " << newCharge->node->siteenergy << " within the interval [" << energypercarrier-0*deltaE << ", " << energypercarrier+2*deltaE << "]" << endl;
                newCharge->node->occupied = 1;
                cout << "starting position for charge " << i + 1 << ": segment " << newCharge->node->id + 1 << endl;
                carrier.push_back(newCharge);
            }

            unsigned insertioncount = 0;
            unsigned long step=0;
            double simtime=0.0;

            vector<int> forbiddennodes;
            vector<int> forbiddendests;
            
            time_t now = time(0);
            tm* localtm = localtime(&now);
            cout << "Run started at " << asctime(localtm) << endl;
 
            double avlifetime=0.0;
            double meanfreepath=0.0;
            vec difflength=vec(0,0,0);
            unsigned steps=0;
            double avgenergy=carrier[0]->node->siteenergy;
            int     carrieridold=carrier[0]->id;
            
            while (insertioncount < insertions) {
                if ((time(NULL) - realtime_start) > _maxrealtime * 60. * 60.) {
                    cout << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime * 60 * 60 + 0.5) << " seconds) has been reached. Stopping here." << endl << endl;
                    break;
                }

                step += 1;
                double cumulated_rate = 0;

                for (unsigned int i = 0; i < carrier.size(); i++) {
                    cumulated_rate += carrier[i]->node->getEscapeRate();
                }
                if (cumulated_rate == 0) { // this should not happen: no possible jumps defined for a node
                    throw runtime_error("ERROR in kmclifetime: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
                }
                // go forward in time
                double dt=Promotetime( cumulated_rate,RandomVariable);
                
                if(_do_carrierenergy){
                    bool print=false;
                    if (carrier[0]->id>carrieridold){
                        avgenergy=carrier[0]->node->siteenergy;
                        print=true;
                        carrieridold=carrier[0]->id;
                    }
                    else if(step%_outputsteps==0){
                        avgenergy=_alpha*carrier[0]->node->siteenergy+(1-_alpha)*avgenergy;
                        print=true;
                    }
                    if(print){
                        energyfile << simtime<<"\t"<<steps <<"\t"<<carrier[0]->id<<"\t"<<avgenergy<<endl;                  
                    }
                }
                simtime += dt;
                steps++;
                for (unsigned int i = 0; i < carrier.size(); i++) {
                    carrier[i]->updateLifetime(dt);
                    carrier[i]->updateSteps(1);
                    carrier[i]->node->occupationtime += dt;
                }
                
                ResetForbiddenlist(forbiddennodes);
                bool secondlevel=true;
                while (secondlevel){

                    // determine which carrier will escape
                    GNode* oldnode=NULL;
                    GNode* newnode=NULL;
                    Chargecarrier* affectedcarrier;

                    double u = 1 - RandomVariable->rand_uniform();
                    for (unsigned int i = 0; i < numberofcharges; i++) {
                        u -= carrier[i]->node->getEscapeRate() / cumulated_rate;
                       
                        if (u <= 0 || i==numberofcharges-1) {
                            oldnode = carrier[i]->node;
                            affectedcarrier = carrier[i];
                            break;}  
                       
                    }
                    vec *dr= NULL;
                   
                    if (CheckForbidden(oldnode->id, forbiddennodes)) {
                        continue;
                    }

                    // determine where it will jump to
                    ResetForbiddenlist(forbiddendests);
                    
                    while (true) {
                        // LEVEL 2
 
                        newnode = NULL;
                        GLink* event=ChooseHoppingDest(oldnode, RandomVariable);
                       
                        if (event->decayevent){
                            oldnode->occupied = 0;
                            avlifetime+=affectedcarrier->getLifetime();
                            meanfreepath+=abs(affectedcarrier->dr_travelled);
                            difflength+=elementwiseproduct(affectedcarrier->dr_travelled,affectedcarrier->dr_travelled);
                            traj << simtime<<"\t"<<insertioncount<< "\t"<< affectedcarrier->id<<"\t"<< affectedcarrier->getLifetime()<<"\t"<<affectedcarrier->getSteps()<<"\t"<< (oldnode->id)+1<<"\t"<<affectedcarrier->dr_travelled.getX()<<"\t"<<affectedcarrier->dr_travelled.getY()<<"\t"<<affectedcarrier->dr_travelled.getZ()<<endl;
                            if(insertions<1500 ||insertioncount% (insertions/1000)==0 || insertioncount<0.001*insertions){
                            std::cout << "\rInsertion " << insertioncount<<" of "<<insertions;
                            std::cout << std::flush;
                            }
                            do {
                                affectedcarrier->node = node[RandomVariable->rand_uniform_int(node.size())];     
                            }
                            while (affectedcarrier->node->occupied == 1 || affectedcarrier->node->injectable != 1 );
                            affectedcarrier->node->occupied=1;
                            affectedcarrier->resetCarrier();
                            insertioncount++;
                            affectedcarrier->id=numberofcharges-1+insertioncount;
                            secondlevel=false;
                            break;
                                }
                        else{
                        newnode = node[event->destination];
                        dr = &(event->dr);
                        }
                    
                        // check after the event if this was allowed
                        if (CheckForbidden(newnode->id, forbiddendests)) {
                            continue;
                        }

                        // if the new segment is unoccupied: jump; if not: add to forbidden list and choose new hopping destination
                        if (newnode->occupied == 1) {
                            if (CheckSurrounded(oldnode, forbiddendests)) {     
                                AddtoForbiddenlist(oldnode->id, forbiddennodes);
                                break; // select new escape node (ends level 2 but without setting level1step to 1)
                            }
                            AddtoForbiddenlist(newnode->id, forbiddendests);
                            continue; // select new destination
                        } else {
                            newnode->occupied = 1;
                            oldnode->occupied = 0;
                            affectedcarrier->node = newnode;
                            affectedcarrier->dr_travelled += *(dr);
                            secondlevel=false;

                            break; // this ends LEVEL 2 , so that the time is updated and the next MC step started
                        }

                       
                        // END LEVEL 2
                    }
                    // END LEVEL 1
                }
            }
            
            
            cout<<endl;
            cout << "Total runtime:\t\t\t\t\t"<< simtime << " s"<< endl;
            cout << "Total KMC steps:\t\t\t\t"<< steps << endl;
            cout << "Average lifetime:\t\t\t\t"<<avlifetime/insertions<< " s"<<endl;
            cout << "Mean freepath\t l=<|r_x-r_o|> :\t\t"<<(meanfreepath/insertions)<< " nm"<<endl;
            cout << "Average diffusionlength\t d=sqrt(<(r_x-r_o)^2>)\t"<<sqrt(abs(difflength)/insertions)<< " nm"<<endl;
            cout<<endl;
            
     for(unsigned int j=0; j<node.size(); j++)
    {   
        occP[j] = node[j]->occupationtime / simtime;
    }
            
            return occP;
        }

void KMCLifetime::WriteOcc(vector<double> occP, vector<GNode*> node)
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

        bool KMCLifetime::EvaluateFrame() {
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "      KMCLIFETIME started" << std::endl;
            std::cout << "-----------------------------------" << std::endl << std::endl;

            // Initialise random number generator
            if (votca::tools::globals::verbose) {
                cout << endl << "Initialising random number generator" << endl;
            }
            srand(_seed); // srand expects any integer in order to initialise the random number generator
            votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
            RandomVariable->init(rand(), rand(), rand(), rand());

            vector<GNode*> node;
            node = LoadGraph();

            InitialRates(node);

            vector<double> occP(node.size(), 0.);

            occP=RunVSSM(node, _insertions, _numberofcharges, RandomVariable);
            WriteOcc(occP, node);
            time_t now = time(0);
            tm* localtm = localtime(&now);
            std::cout << "      KMCLIFETIME finished at:" <<asctime(localtm) <<  std::endl;
            
        

            return true;
        }
    }
}
