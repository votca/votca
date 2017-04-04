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



#include <votca/tools/tokenizer.h>
#include <votca/xtp/kmccalculator.h>

#include <votca/tools/constants.h>
namespace votca { namespace xtp {

  
   


class KMCMultiple : public KMCCalculator 
{
public:
    KMCMultiple() {};
   ~KMCMultiple() {};
   std::string Identify() { return "kmcmultiple"; }
    void Initialize(tools::Property *options);
     bool EvaluateFrame(ctp::Topology *top);



private:
        
            void  RunVSSM(ctp::Topology *top);
            int _maxsteps;
            double _runtime;
            double _outputtime;
            std::string _trajectoryfile;
            std::string _timefile;
            double _maxrealtime;
            bool    _tof;
            double _toflength;
            tools::vec _tofdirection;
};




void KMCMultiple::Initialize(tools::Property *options){
    std::string key = "options." + Identify();

    _runtime=options->ifExistsReturnElseThrowRuntimeError<double>(key+".runtime");
    _seed=options->ifExistsReturnElseThrowRuntimeError<int>(key+"seed");
    _numberofcharges=options->ifExistsReturnElseThrowRuntimeError<int>(key+".numberofcharges");
    _injection_name=options->ifExistsReturnElseThrowRuntimeError<int>(key+".injectionpattern");
  

        _maxrealtime=options->ifExistsReturnElseReturnDefault<double>(key+".maxrealtime",1E10);
        _trajectoryfile=options->ifExistsReturnElseReturnDefault<std::string>(key+".trajectoryfile","trajectory.csv");
        _temperature=options->ifExistsReturnElseReturnDefault<double>(key+".temperature",300);
        _rates=options->ifExistsReturnElseReturnDefault<std::string>(key+".rates","statefile");
        
     
        _injectionmethod = options->ifExistsReturnElseReturnDefault<std::string>(key+".injectionmethod","random");
	
        if (_injectionmethod != "random" && _injectionmethod != "equilibrated")
        {
	    cout << "WARNING in kmcmultiple: Unknown injection method. It will be set to random injection." << endl;
            _injectionmethod = "random";
            
        }
         _field = options->ifExistsReturnElseReturnDefault<tools::vec>(key+".field",tools::vec(0,0,0));
        
      
	_outputtime = options->ifExistsReturnElseReturnDefault<double>(key+".outputtime",0);
        _timefile = options->ifExistsReturnElseReturnDefault<std::string>(key+".timefile","timedependence.csv");
	
        std::string carriertype=options->ifExistsReturnElseReturnDefault<std::string>(key+".carriertype","e");
        _carriertype=StringtoCarriertype(carriertype);
     
        


        _tof = options->ifExistsReturnElseReturnDefault<bool>(key+".tof",false);
        _tofdirection = options->ifExistsReturnElseReturnDefault<tools::vec>(key+".tofdirection",tools::vec(0,1,0));
         _tofdirection=_tofdirection.normalize();
        _toflength =options->ifExistsReturnElseReturnDefault<double>(key+".toflength",50);

        
        if(_tof){
            cout << "Time of flight experiment: ON" << endl;
            cout << "direction for TOF condition: " << _tofdirection << endl;
            cout << "Sample length for TOF experiment: " << _toflength<< " nm" << endl;
        }
        else{
            cout << "Time of flight experiment: OFF (bulk mode)" << endl;
        }

        return;
}      
        


void KMCMultiple::RunVSSM(ctp::Topology *top)
{

    int realtime_start = time(NULL);
    cout << endl << "Algorithm: VSSM for Multiple Charges" << endl;
    cout << "number of charges: " << _numberofcharges << endl;
    cout << "number of nodes: " << _nodes.size() << endl;
    string stopcondition;
    unsigned long maxsteps=0;
    int diffusionsteps = 0;
    unsigned diffusion_stepsize = 10000;
    tools::matrix avgdiffusiontensor;
    avgdiffusiontensor.ZeroMatrix();
    double accumulatedenergy = 0;
    int tdpendencesteps = 1;
    double currentmobility_laststep = 7E77;
    int currentkmcstep_laststep = 0;
    double relaxationtime=0.0;
    double relaxationlength=0.0;
    
    int mobilitycheck = 0;
    if(_runtime > 100)
    {
        stopcondition = "steps";
        maxsteps = _runtime;
        cout << "stop condition: " << maxsteps << " steps." << endl;
    }
    else
    {
        stopcondition = "runtime";
        cout << "stop condition: " << _runtime << " seconds runtime." << endl;
        cout << "(If you specify runtimes larger than 100 kmcmultiple assumes that you are specifying the number of steps.)" << endl;
    }
    
    if(_numberofcharges > _nodes.size())
    {
        throw runtime_error("ERROR in kmcmultiple: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
    }
    
    fstream traj;
   
    cout << "Writing trajectory to " << _trajectoryfile << "." << endl; 
    traj.open (_trajectoryfile.c_str(), fstream::out);
    fstream tfile;
    cout << "Writing time dependence of energy and mobility to " << _timefile << "." << endl; 
    tfile.open (_timefile.c_str(), fstream::out);
    if(_outputtime != 0)
    {   
        traj << "'time[s]'\t";
        for(unsigned int i=0; i<_numberofcharges; i++)
        {
            traj << "'carrier" << i+1 << "_x'\t";    
            traj << "'carrier" << i+1 << "_y'\t";    
            traj << "'carrier" << i+1 << "_z";    
            if(i<_numberofcharges-1)
            {
                traj << "'\t";
            }
        }
        traj << endl;
        
        tfile << "'time[s]'\t'energy_per_carrier[eV]\t'mobility[m**2/Vs]'\t'distance_fielddirection[m]'\t'distance_absolute[m]'\t'diffusion_coefficient[m**2]'" << endl;
        
    }
    double outputfrequency = _runtime/100;
    double outputstepfrequency = _maxsteps/100;
    
    
    double absolute_field = tools::abs(_field);

    // Injection
    cout << endl << "injection method: " << _injectionmethod << endl;
    double deltaE = 0;
    double energypercarrier=0;
    double totalenergy=0;
    if(_injectionmethod == "equilibrated"){
        vector< double > energy;
        for(unsigned int i=0; i<_nodes.size(); i++)
        {
            energy.push_back(_nodes[i]->siteenergy);
        }
        std::sort (energy.begin(), energy.end());
        double fermienergy = 0.5*(energy[_numberofcharges-1]+energy[_numberofcharges]);
                cout << "Energy range in morphology data: [" << energy[0] << ", " << energy[energy.size()-1] << "] eV." << endl;
        cout << "first guess Fermi energy: " << fermienergy << " eV."<< endl;

        cout << "improving guess for fermi energy ";
        
        double probsum;
        double fstepwidth = std::abs(0.01 * fermienergy);
        probsum = 0;
        while(std::abs(probsum - _numberofcharges) > 0.1)
        {
            cout << ".";
            totalenergy = 0;
            probsum = 0;
            for(unsigned int i=0; i<energy.size(); i++)
            {
                totalenergy += energy[i] * 1/(exp((energy[i]-fermienergy)/tools::conv::kB/_temperature)+1);
                probsum += 1/(exp((energy[i]-fermienergy)/tools::conv::kB/_temperature)+1);
                
            }
            if(std::abs(probsum) > _numberofcharges)
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
        while(inside < _numberofcharges){
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
    vector< Chargecarrier* > carriers;
    vector<tools::vec> startposition(_numberofcharges,tools::vec(0,0,0));
    cout << "looking for injectable nodes..." << endl;
    for (unsigned int i = 0; i < _numberofcharges; i++){
        Chargecarrier *newCharge = new Chargecarrier;      
        newCharge->id = i;
        newCharge->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];
        int ininterval = 1;
        if (_injectionmethod == "equilibrated") {ininterval = 0;}
        while(newCharge->node->occupied == 1 || newCharge->node->injectable != 1 || ininterval != 1)
        {   // maybe already occupied? or maybe not injectable?
            newCharge->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];
            if(_injectionmethod == "equilibrated")
            { // check if charge is in the acceptance interval
                if(newCharge->node->siteenergy >= energypercarrier-0*deltaE && 
                   newCharge->node->siteenergy <= energypercarrier+2*deltaE && newCharge->node->occupied == 0 && newCharge->node->injectable == 1)
                {
                    ininterval = 1;
                }
                
            }
        }
        // cout << "selected segment " << newCharge->node->id+1 << " which has energy " << newCharge->node->siteenergy << " within the interval [" << energypercarrier-0*deltaE << ", " << energypercarrier+2*deltaE << "]" << endl;
        newCharge->node->occupied = 1;
        newCharge->resetCarrier();
        startposition[i] = newCharge->node->position;
        cout << "starting position for charge " << i+1 << ": segment " << newCharge->node->id+1 << endl;
        carriers.push_back(newCharge);
    }
    
    if(_injectionmethod == "equilibrated")
    {
        cout << "repositioning charges to obtain an equilibrium configuration..." << endl;
        double realisedenergy = 0;
        //while(std::abs(realisedenergy - totalenergy) > 0.01*(std::abs(realisedenergy)+std::abs(totalenergy)))
        //{
            realisedenergy = 0;
            for (unsigned int j = 0; j < _numberofcharges; j++)
            {
                realisedenergy += carriers[j]->node->siteenergy;
            }

            for(unsigned int i=0; i<_numberofcharges; i++)
            {
                for(unsigned int k=0; k<_nodes.size(); k++)
                {
                    if(std::abs(realisedenergy-carriers[i]->node->siteenergy+_nodes[k]->siteenergy - totalenergy) < std::abs(realisedenergy - totalenergy) && _nodes[k]->occupied == 0)
                    {    //move charge
                         carriers[i]->node->occupied = 0;
                         realisedenergy = realisedenergy-carriers[i]->node->siteenergy+_nodes[k]->siteenergy;
                         carriers[i]->node = _nodes[k];
                         _nodes[k]->occupied = 1;
                    }        
                }
                
            }
        //}    
        cout << "realised energy: " << realisedenergy/_numberofcharges << " eV (aimed for " << energypercarrier << " eV)"<<  endl;
        for(unsigned int i=0; i<_numberofcharges; i++)
        {
            startposition[i] = carriers[i]->node->position;
            cout << "starting position for charge " << i+1 << ": segment " << carriers[i]->node->id+1 << endl;
        }
    }
    

    double simtime = 0.;
    unsigned long step = 0;
    double nextoutput = outputfrequency;
    unsigned long nextstepoutput = outputstepfrequency;
    double nexttrajoutput = _outputtime;
    unsigned nextdiffstep = diffusion_stepsize;
    
    vector<int> forbiddennodes;
    vector<int> forbiddendests;
    
    bool timeout = false;
    while(((stopcondition == "runtime" && simtime < _runtime) || (stopcondition == "steps" && step < maxsteps)) && timeout == false)
    {
        if((time(NULL) - realtime_start) > _maxrealtime*60.*60.)
        {
            cout  << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime*60*60 +0.5) <<" seconds) has been reached. Stopping here." << endl << endl;
            break;
        }
        double cumulated_rate = 0;
       
        for(unsigned int i=0; i<_numberofcharges; i++)
        {
            cumulated_rate += carriers[i]->node->getEscapeRate();
        }
        if(cumulated_rate == 0)
        {   // this should not happen: no possible jumps defined for a node
            throw runtime_error("ERROR in kmcmultiple: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
        }
        
        double dt =Promotetime(cumulated_rate);
        
        simtime += dt;
        if(tools::globals::verbose) {cout << "simtime += " << dt << endl << endl;}
        step += 1;
        
        for(unsigned int i=0; i<_numberofcharges; i++)
        {
            carriers[i]->node->occupationtime += dt;
        }

        
        ResetForbiddenlist(forbiddennodes);
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
            
            double u = 1 - _RandomVariable->rand_uniform();
            for(unsigned int i=0; i<_numberofcharges; i++)
            {
                u -= carriers[i]->node->getEscapeRate()/cumulated_rate;
                if(u <= 0)
                {
                    do_oldnode = carriers[i]->node;
                    do_affectedcarrier = carriers[i];
                    break;
                }
               do_oldnode = carriers[i]->node;
               do_affectedcarrier = carriers[i];
            }
                
            //double maxprob = 0.;
            //double newprob = 0.;
            tools::vec dr;
            if(tools::globals::verbose) {cout << "Charge number " << do_affectedcarrier->id+1 << " which is sitting on segment " << do_oldnode->id+1 << " will escape!" << endl ;}
            if(CheckForbidden(do_oldnode->id, forbiddennodes) == 1) {continue;}
            
            // determine where it will jump to
            ResetForbiddenlist(forbiddendests);
            while(timeout == false)
            {
            // LEVEL 2
                if(tools::globals::verbose) {cout << "There are " << do_oldnode->events.size() << " possible jumps for this charge:"; }
                
                if((time(NULL) - realtime_start) > _maxrealtime*60.*60.)
                {
                    cout  << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime*60*60 +0.5) <<" seconds) has been reached while searching destination node (level 2). Stopping here." << endl << endl;
                    timeout = true;
                    break;
                }


                do_newnode = NULL;
                u = 1 - _RandomVariable->rand_uniform();
                for(unsigned int j=0; j<do_oldnode->events.size(); j++){
                    if(tools::globals::verbose) { cout << " " << do_oldnode->events[j].destination+1 ; }
                    u -= do_oldnode->events[j].rate/do_oldnode->getEscapeRate();
                    if(u <= 0){
                        do_newnode = _nodes[do_oldnode->events[j].destination];
                        dr = do_oldnode->events[j].dr;
                        break;
                    }
                    do_newnode = _nodes[do_oldnode->events[j].destination];
                    dr = do_oldnode->events[j].dr;
                }

                if(do_newnode == NULL){
                    if(tools::globals::verbose) {
                        cout << endl << "Node " << do_oldnode->id+1  << " is SURROUNDED by forbidden destinations and zero rates. "
                                "Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl; 
                    }
                    AddtoForbiddenlist(do_oldnode->id, forbiddennodes);
                    break; // select new escape node (ends level 2 but without setting level1step to 1)
                }
                if(tools::globals::verbose) {cout << endl << "Selected jump: " << do_newnode->id+1 << endl; }
                
                // check after the event if this was allowed
                if(CheckForbidden(do_newnode->id, forbiddendests) == 1){
                    if(tools::globals::verbose) {cout << "Node " << do_newnode->id+1  << " is FORBIDDEN. Now selection new hopping destination." << endl; }
                    continue;
                }

                // if the new segment is unoccupied: jump; if not: add to forbidden list and choose new hopping destination
                if(do_newnode->occupied == 1){
                    if(CheckSurrounded(do_oldnode, forbiddendests)){
                        if(tools::globals::verbose) {
                            cout << "Node " << do_oldnode->id+1  << " is SURROUNDED by forbidden destinations. "
                                    "Adding it to the list of forbidden nodes. After that: selection of a new escape node." << endl; 
                        }
                        AddtoForbiddenlist(do_oldnode->id, forbiddennodes);
                        break; // select new escape node (ends level 2 but without setting level1step to 1)
                    }
                    if(tools::globals::verbose) {cout << "Selected segment: " << do_newnode->id+1 << " is already OCCUPIED. Added to forbidden list." << endl << endl;}
                    AddtoForbiddenlist(do_newnode->id, forbiddendests);
                    if(tools::globals::verbose) {cout << "Now choosing different hopping destination." << endl; }
                    continue; // select new destination
                }
                else{
                    do_newnode->occupied = 1;
                    do_oldnode->occupied = 0;
                    do_affectedcarrier->node = do_newnode;
                    do_affectedcarrier->dr_travelled += dr;
                    level1step = 1;
                    if(tools::globals::verbose) {cout << "Charge has jumped to segment: " << do_newnode->id+1 << "." << endl;}
                    
                    if (_tof){   // time-of-flight experiment
                        
                        // keep track of travel in TOF direction
                        
                        
                        
                        do_affectedcarrier->dr_travelled += dr;
                        double dr_tof = _tofdirection*do_affectedcarrier->dr_travelled;
                       
                        
                        // boundary crossing
                        if(dr_tof > _toflength){
                            cout << "\nTOF length has been crossed in "<<_tofdirection<<" direction." << endl;
                            cout << "LAST POSITION: " << (do_affectedcarrier->dr_travelled-dr) << ", NEW POSITION: " << do_affectedcarrier->dr_travelled << endl ; 

                            // now re-inject the carrier to a new random position
                            GNode *temp_node = new GNode;
                            temp_node = do_affectedcarrier->node;
                            do_affectedcarrier->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];
                            while(do_affectedcarrier->node->occupied == 1 || do_affectedcarrier->node->injectable != 1)
                            {   // maybe already occupied? or maybe not injectable?
                                do_affectedcarrier->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];
                            }
                            do_affectedcarrier->node->occupied = 1;
                            temp_node->occupied = 0;
                            cout << "carrier " << do_affectedcarrier->id << " has been instanteneously moved from node " << temp_node->id+1 << " to node " << do_affectedcarrier->node->id+1 << endl << endl; 
                            
                            // reset TOF length counter
                            do_affectedcarrier->resetCarrier();
                        }
                        
                    }
                    
                    break; // this ends LEVEL 2 , so that the time is updated and the next MC step started
                }

                if(tools::globals::verbose) {cout << "." << endl;}
            // END LEVEL 2
            }
        // END LEVEL 1
        }    
        
        
        if(step > nextdiffstep)       
        {
            nextdiffstep += diffusion_stepsize;
            for(unsigned int i=0; i<_numberofcharges; i++) 
            {
                diffusionsteps  += 1;
                avgdiffusiontensor += (carriers[i]->dr_travelled)|(carriers[i]->dr_travelled);
            }
        }
        
        if(_outputtime != 0 && simtime > nexttrajoutput)       
        {
            // write to trajectory file
            nexttrajoutput = simtime + _outputtime;
            traj << simtime << "\t";
            for(unsigned int i=0; i<_numberofcharges; i++) 
            {
                traj << startposition[i].getX() + carriers[i]->dr_travelled.getX() << "\t";
                traj << startposition[i].getY() + carriers[i]->dr_travelled.getY() << "\t";
                traj << startposition[i].getZ() + carriers[i]->dr_travelled.getZ();
                if (i<_numberofcharges-1) 
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
            tools::vec dr_travelled_current = tools::vec (0,0,0);
            if(absolute_field != 0)
            {
                tools::vec avgvelocity_current = tools::vec(0,0,0);
                for(unsigned int i=0; i<_numberofcharges; i++)
                {
                    dr_travelled_current += carriers[i]->dr_travelled;
                    accumulatedenergy += carriers[i]->node->siteenergy;
                }
                dr_travelled_current /= _numberofcharges;
                currentenergy = accumulatedenergy /tdpendencesteps/_numberofcharges;
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
                for(unsigned int i=0; i<_numberofcharges; i++)
                {
                   
                    relaxationlength +=abs(carriers[i]->dr_travelled);
                }
                relaxationlength /= _numberofcharges;
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
        
           
        }
        else if(stopcondition == "steps" && step > nextstepoutput)
        {
            nextstepoutput = step + outputstepfrequency;
           
            
          
        }
    }
    
    if(_outputtime != 0)
    {   
        traj.close();
        tfile.close();
    }


    
    vector< ctp::Segment* >& seg = top->Segments();
    for (unsigned i = 0; i < seg.size(); i++) {
            double occupationprobability=_nodes[i]->occupationtime / simtime;
            seg[i]->setOcc(occupationprobability,_carriertype);
        }

    cout << endl << "finished KMC simulation after " << step << " steps." << endl;
    cout << "simulated time " << simtime << " seconds." << endl;
    cout << "runtime: ";
    cout << endl << endl;
    tools::vec dr_travelled = tools::vec (0,0,0);
    tools::vec avgvelocity = tools::vec(0,0,0);
    for(unsigned int i=0; i<_numberofcharges; i++)
    {
        cout << std::scientific << "    charge " << i+1 << ": " << carriers[i]->dr_travelled/simtime << endl;
        dr_travelled += carriers[i]->dr_travelled;
    }
    dr_travelled /= _numberofcharges;
    
    avgvelocity = dr_travelled/simtime; 
    cout << std::scientific << "  Overall average velocity (nm/s): " << avgvelocity << endl;

    cout << endl << "Distances travelled (nm): " << endl;
    for(unsigned int i=0; i<_numberofcharges; i++)
    {
        cout << std::scientific << "    charge " << i+1 << ": " << carriers[i]->dr_travelled << endl;
    }
    
    // calculate mobilities
    double average_mobility = 0;
    if (absolute_field != 0)
    {
        cout << endl << "Mobilities (nm^2/Vs): " << endl;
        for(unsigned int i=0; i<_numberofcharges; i++)
        {
            //tools::vec velocity = carrier[i]->dr_travelled/simtime*1e-9;
            tools::vec velocity = carriers[i]->dr_travelled/simtime;
            //double absolute_velocity = sqrt(velocity.x()*velocity.x() + velocity.y()*velocity.y() + velocity.z()*velocity.z());
            //cout << std::scientific << "    charge " << i+1 << ": mu=" << absolute_velocity/absolute_field*1E4 << endl;
            cout << std::scientific << "    charge " << i+1 << ": mu=" << (velocity*_field)/absolute_field/absolute_field << endl;
            average_mobility += (velocity*_field) /absolute_field/absolute_field;
        }
        average_mobility /= _numberofcharges;
        cout << std::scientific << "  Overall average mobility in field direction <mu>=" << average_mobility << " m^2/Vs (= " << average_mobility*1E4 << "cm^2/Vs)" << endl;
      }
    cout << endl;
    
    // calculate diffusion tensor
    avgdiffusiontensor /= (diffusionsteps*2*simtime*_numberofcharges);
    cout<<endl<<"Diffusion tensor averaged over all carriers (nm^2/s):" << endl << avgdiffusiontensor << endl;

    tools::matrix::eigensystem_t diff_tensor_eigensystem;
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
       average_mobility = std::abs(avgD / tools::conv::kB / _temperature);
       cout << std::scientific << "  Overall average mobility <mu>=" << average_mobility << " nm^2/Vs (= "  << endl;
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


    
    return;
}




bool KMCMultiple::EvaluateFrame(ctp::Topology *top){
    std::cout << "-----------------------------------" << std::endl;      
    std::cout << "      KMC FOR MULTIPLE CHARGES" << std::endl;
    std::cout << "-----------------------------------" << std::endl << std::endl;      
 
    // Initialise random number generator
    if(tools::globals::verbose) { cout << endl << "Initialising random number generator" << endl; }
    srand(_seed); // srand expects any integer in order to initialise the random number generator
    _RandomVariable = new tools::Random2();
    _RandomVariable->init(rand(), rand(), rand(), rand());
    
    LoadGraph(top);
    
        
        if(_rates == "calculate")
        {
            cout << "Calculating rates (i.e. rates from state file are not used)." << endl;
            InitialRates();
        }
        else
        {
            cout << "Using rates from state file." << endl;
        }
    

    RunVSSM(top);

    return true;
}

}}


#endif	/* __VOTCA_KMC_MULTIPLE_H */
