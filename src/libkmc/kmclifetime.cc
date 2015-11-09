/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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


using namespace std;

namespace votca {
    namespace kmc {

        const static double kB = 8.617332478E-5; // eV/K
        const static double hbar = 6.5821192815E-16; // eV*s
        const static double eps0 = 8.85418781762E-12 / 1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
        const static double epsr = 3.0; // relative material permittivity
        const static double Pi = 3.14159265358979323846;

        void KMCLifetime::Initialize(const char *filename, Property *options, const char *outputfile) {
            if (options->exists("options.kmclifetime.insertions")) {
                _insertions = options->get("options.kmclifetime.insertions").as<double>();
            } else {
                throw runtime_error("Error in kmclifetime: total number of insertions not provided");
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
                _lifetimefile = options->get("options.kmclifetime.lifetimes").as<string>();
            } else {
                throw runtime_error("Error in kmclifetime: injection pattern is not provided");
            }

            if (options->exists("options.kmclifetime.trajectoryfile")) {
                _trajectoryfile = options->get("options.kmclifetime.trajectoryfile").as<string>();
            } else {
                _trajectoryfile = "trajectory.csv";
            }
            if (_trajectoryfile == "") {
                _trajectoryfile = "trajectory.csv";
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
                cout << "WARNING in kmclifetime: You did not specify a temperature. If no explicit Coulomb interaction is used, this is not a problem, as the rates are read from the state file and temperature is not needed explicitly in the KMC simulation. Otherwise a default value of 300 K is used." << endl;
                _temperature = 0;
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

        void KMCLifetime::progressbar(double fraction) {
            int totalbars = 50;
            std::cout << "\r";
            for (double bars = 0; bars<double(totalbars); bars++) {
                if (bars <= fraction * double(totalbars)) {
                    std::cout << "|";
                } else {
                    std::cout << "-";
                }
            }
            std::cout << "  " << int(fraction * 1000) / 10. << " %   ";
            std::cout << std::flush;
            if (fraction * 100 == 100) {
                std::cout << std::endl;
            }
        }

        
        
        
        vector<GNode*> KMCLifetime::LoadGraph() {
            vector<GNode*> node;

            GNode *newNode = new GNode();
            
            node.push_back(newNode);
            
            
            // Load nodes
            votca::tools::Database db;
            db.Open(_filename);
            if (votca::tools::globals::verbose) {
                cout << "LOADING GRAPH" << endl << "database file: " << _filename << endl;
            }
            votca::tools::Statement *stmt;
            if (_carriertype == "h" || _carriertype == "e") {
                stmt = db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnCnN" + _carriertype + ", UcNcC" + _carriertype + ",eAnion,eNeutral,eCation,UcCnN" + _carriertype + " FROM segments;");
            } else if (_carriertype == "s" || _carriertype == "t") {
                stmt = db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnXnN" + _carriertype + ", UxNxX" + _carriertype + ",eSinglet,eNeutral,eTriplet,UxXnN" + _carriertype + " FROM segments;");
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
                double eNeutral = stmt->Column<double>(8);
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
            if (votca::tools::globals::verbose) {
                cout << "segments: " << node.size() << endl;
            }

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

            //blabla

            // Calculate initial escape rates !!!THIS SHOULD BE MOVED SO THAT IT'S NOT DONE TWICE IN CASE OF COULOMB INTERACTION!!!
            for (unsigned int i = 0; i < node.size(); i++) {
                node[i]->InitEscapeRate();
            }
            return node;
        }

        void KMCLifetime::ResetForbidden(vector<int> &forbiddenid) {
            forbiddenid.clear();
        }

        void KMCLifetime::AddForbidden(int id, vector<int> &forbiddenid) {
            forbiddenid.push_back(id);
        }

        bool KMCLifetime::Forbidden(int id, vector<int> forbiddenlist) {
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
/*
        bool KMCLifetime::Surrounded(GNode* node, vector<int> forbiddendests) {
            bool surrounded = true;
            for (unsigned int i = 0; i < node->event.size(); i++) {
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
*/
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
            int n = sprintf(buffer, "%d:%02d:%02d", hours, minutes, seconds);
            printf("%s", buffer, n);
        }

        void KMCLifetime::InitialRates(vector<GNode*> node) {
            cout << endl << "Calculating initial Marcus rates." << endl;
            cout << "    Temperature T = " << _temperature << " K." << endl;
           
            cout << "    carriertype: " << _carriertype << endl;
            int numberofsites = node.size();
            cout << "    Rates for " << numberofsites << " sites are computed." << endl;
            double maxreldiff = 0;
            int totalnumberofrates = 0;
            for (unsigned int i = 0; i < numberofsites; i++) {
                int numberofneighbours = node[i]->event.size();
                for (unsigned int j = 0; j < numberofneighbours; j++) {
                    double dX = node[i]->event[j].dr.x();
                    double dY = node[i]->event[j].dr.y();
                    double dZ = node[i]->event[j].dr.z();
                    

                    double destindex = node[i]->event[j].destination;
                    double reorg = node[i]->reorg_intorig + node[destindex]->reorg_intdest + node[i]->event[j].reorg_out;

                    double dG = node[destindex]->siteenergy - node[i]->siteenergy;

                    double J2 = node[i]->event[j].Jeff2;

                    double rate = 2 * Pi / hbar * J2 / sqrt(4 * Pi * reorg * kB * _temperature) * exp(-(dG + reorg)*(dG + reorg) / (4 * reorg * kB * _temperature));

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
        }

        void KMCLifetime::RunVSSM(vector<GNode*> node, double runtime, unsigned int numberofcharges, votca::tools::Random2 * RandomVariable) {

            int realtime_start = time(NULL);
            cout << endl << "Algorithm: VSSM for Multiple Charges" << endl;
            cout << "number of charges: " << numberofcharges << endl;
            cout << "number of nodes: " << node.size() << endl;
            

            if (numberofcharges > node.size()) {
                throw runtime_error("ERROR in kmclifetime: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
            }

            fstream traj;
            char trajfile[100];
            strcpy(trajfile, _trajectoryfile.c_str());
            cout << "Writing trajectory to " << trajfile << "." << endl;
            traj.open(trajfile, fstream::out);
            if (_outputtime != 0) {
                traj << "'time[s]'\t";
                for (unsigned int i = 0; i < numberofcharges; i++) {
                    traj << "'carrier" << i + 1 << "_x'\t";
                    traj << "'carrier" << i + 1 << "_y'\t";
                    traj << "'carrier" << i + 1 << "_z";
                    if (i < numberofcharges - 1) {
                        traj << "'\t";
                    }
                }
                traj << endl;

              

            }
            
            vector<double> occP(node.size(), 0.);

            

            // Injection
            cout << endl << "injection method: " << _injectionmethod << endl;
            
            vector< Chargecarrier* > carrier;
            vector<vec> startposition(numberofcharges, vec(0, 0, 0));
            cout << "looking for injectable nodes..." << endl;
            for (unsigned int i = 0; i < numberofcharges; i++) {
                Chargecarrier *newCharge = new Chargecarrier;
                newCharge->id = i;
                newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];
                int ininterval = 1;
               
                while (newCharge->node->occupied == 1 || newCharge->node->injectable != 1 || ininterval != 1) { // maybe already occupied? or maybe not injectable?
                    newCharge->node = node[RandomVariable->rand_uniform_int(node.size())];  
                }
                // cout << "selected segment " << newCharge->node->id+1 << " which has energy " << newCharge->node->siteenergy << " within the interval [" << energypercarrier-0*deltaE << ", " << energypercarrier+2*deltaE << "]" << endl;
                newCharge->node->occupied = 1;
               
     
                
                cout << "starting position for charge " << i + 1 << ": segment " << newCharge->node->id + 1 << endl;
                carrier.push_back(newCharge);
            }

         

            int realtime_start = time(NULL);

            int insertions = 0;
            unsigned long step=0;
            double simtime=0.0;



            progressbar(0.);
            vector<int> forbiddennodes;
            vector<int> forbiddendests;

            while (insertions < _insertions) {
                if ((time(NULL) - realtime_start) > _maxrealtime * 60. * 60.) {
                    cout << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime * 60 * 60 + 0.5) << " seconds) has been reached. Stopping here." << endl << endl;
                    break;
                }
                
                double cumulated_rate = 0;

                for (unsigned int i = 0; i < carrier.size(); i++) {
                    cumulated_rate += carrier[i]->node->getEscapeRate();
                }
                if (cumulated_rate == 0) { // this should not happen: no possible jumps defined for a node
                    throw runtime_error("ERROR in kmclifetime: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
                }
                // go forward in time
                double dt = 0;
                double rand_u = 1 - RandomVariable->rand_uniform();
                while (rand_u == 0) {
                    cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                    rand_u = 1 - RandomVariable->rand_uniform();
                }
                dt = -1 / cumulated_rate * log(rand_u);
                simtime += dt;
                
                for (unsigned int i = 0; i < carrier.size(); i++) {
                    cumulated_rate += carrier[i]->updateLifetime(dt);
                }
                if (votca::tools::globals::verbose) {
                    cout << "simtime += " << dt << endl << endl;
                }
                step += 1;

                for (unsigned int i = 0; i < numberofcharges; i++) {
                    carrier[i]->node->occupationtime += dt;
                }


                ResetForbidden(forbiddennodes);
                int level1step = 0;
                while (level1step == 0 )
                    // LEVEL 1
                {



                    // determine which electron will escape
                    GNode* do_oldnode;
                    GNode* do_newnode;
                    Chargecarrier* do_affectedcarrier;

                    double u = 1 - RandomVariable->rand_uniform();
                    for (unsigned int i = 0; i < numberofcharges; i++) {
                        u -= carrier[i]->node->EscapeRate() / cumulated_rate;
                        if (u <= 0) {
                            do_oldnode = carrier[i]->node;
                            do_affectedcarrier = carrier[i];
                            break;
                        }
                        do_oldnode = carrier[i]->node;
                        do_affectedcarrier = carrier[i];
                    }

                  
                    vec dr;
                    if (votca::tools::globals::verbose) {
                        cout << "Charge number " << do_affectedcarrier->id + 1 << " which is sitting on segment " << do_oldnode->id + 1 << " will escape!" << endl;
                    }
                    if (Forbidden(do_oldnode->id, forbiddennodes) == 1) {
                        continue;
                    }

                    // determine where it will jump to
                    ResetForbidden(forbiddendests);
                    while (true) {
                        // LEVEL 2
                        if (votca::tools::globals::verbose) {
                            cout << "There are " << do_oldnode->event.size() << " possible jumps for this charge:";
                        }

                       


                        do_newnode = NULL;
                        u = 1 - RandomVariable->rand_uniform();
                        for (unsigned int j = 0; j < do_oldnode->event.size(); j++) {
                            if (votca::tools::globals::verbose) {
                                cout << " " << do_oldnode->event[j].destination + 1;
                            }
                            u -= do_oldnode->event[j].rate / do_oldnode->EscapeRate();
                            if (u <= 0) {
                                do_newnode = node[do_oldnode->event[j].destination];
                                dr = do_oldnode->event[j].dr;
                                break;
                            }
                            do_newnode = node[do_oldnode->event[j].destination];
                            dr = do_oldnode->event[j].dr;
                        }

                        if (do_newnode == NULL) {
                            if (votca::tools::globals::verbose) {
                                cout << endl << "Node " << do_oldnode->id + 1 << " is SURROUNDED by forbidden destinations and zero rates. Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl;
                            }
                            AddForbidden(do_oldnode->id, forbiddennodes);
                            break; // select new escape node (ends level 2 but without setting level1step to 1)
                        }
                        
                        if (votca::tools::globals::verbose) {
                            cout << endl << "Selected jump: " << do_newnode->id + 1 << endl;
                        }

                        // check after the event if this was allowed
                        if (Forbidden(do_newnode->id, forbiddendests) == true) {
                            if (votca::tools::globals::verbose) {
                                cout << "Node " << do_newnode->id + 1 << " is FORBIDDEN. Now selection new hopping destination." << endl;
                            }
                            continue;
                        }

                        // if the new segment is unoccupied: jump; if not: add to forbidden list and choose new hopping destination
                        if (do_newnode->occupied == 1) {
                            if (Surrounded(do_oldnode, forbiddendests) == true) {
                                if (votca::tools::globals::verbose) {
                                    cout << "Node " << do_oldnode->id + 1 << " is SURROUNDED by forbidden destinations. Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl;
                                }
                                AddForbidden(do_oldnode->id, forbiddennodes);
                                break; // select new escape node (ends level 2 but without setting level1step to 1)
                            }
                            if (votca::tools::globals::verbose) {
                                cout << "Selected segment: " << do_newnode->id + 1 << " is already OCCUPIED. Added to forbidden list." << endl << endl;
                            }
                            AddForbidden(do_newnode->id, forbiddendests);
                            if (votca::tools::globals::verbose) {
                                cout << "Now choosing different hopping destination." << endl;
                            }
                            continue; // select new destination
                        } else {
                            do_newnode->occupied = 1;
                            do_oldnode->occupied = 0;
                            do_affectedcarrier->node = do_newnode;
                            do_affectedcarrier->dr_travelled += dr;
                            level1step = 1;
                            if (votca::tools::globals::verbose) {
                                cout << "Charge has jumped to segment: " << do_newnode->id + 1 << "." << endl;
                            }



                            break; // this ends LEVEL 2 , so that the time is updated and the next MC step started
                        }

                        if (votca::tools::globals::verbose) {
                            cout << "." << endl;
                        }
                        // END LEVEL 2
                    }
                    // END LEVEL 1
                }
            }
            return occP;
        }



        bool KMCLifetime::EvaluateFrame() {
            std::cout << "-----------------------------------" << std::endl;
            std::cout << "      KMCLIFETIME" << std::endl;
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

            RunVSSM(node, _runtime, _numberofcharges, RandomVariable);

        

            return true;
        }
    }
}
