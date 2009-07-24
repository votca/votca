// 
// File:   csg_nemat.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:35 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <tools/cubicspline.h>
#include <neighbourlist.h>
#include <nblist.h>
#include <beadlist.h>
#include <exclusionlist.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <sstream>
#include "csg_fmatch.h"


//#define _DEBUG

//namespace ub = boost::numeric::ublas;
//using namespace std;

void CGForceMatching::BeginCG(Topology *top, Topology *top_atom)
    {

        int lines_init = 0, colms_init = 0;  // initial size of _A 
        int sfnum; // number of spline functions for a cubicspline. 
        
        double grid_min, grid_max, grid_step; // used for initializing grid for 
                                              // different interactions
                
        int interaction_number = 0;
        NewNeighbour = true;
  //      beadTypes = 1; // used by CGForceMatching::beadType2intType
        beadTypes = _options.get("cg.fmatch.bead_types").as<int>();
        numBondInt = 0;
        ConstrLeastSQ = false;
//        N_frames = 6; // Number of frames in the block
        N_frames = _options.get("cg.fmatch.frames_per_block").as<int>();
        BlockNum = 0;
        cutoff = _options.get("cg.fmatch.cutoff").as<double>();

        // set counters to zero value:
        line_cntr = col_cntr = 0;
        
        // initializing bonded interactions
        for (list<Property*>::iterator iter = _bonded.begin();
               iter != _bonded.end(); ++iter) {
            
            SplineInfo *i = new SplineInfo;
            i->splineName = (*iter)->get("name").value();
            grid_min = (*iter)->get("fmatch.min").as<double>();
            grid_max = (*iter)->get("fmatch.max").as<double>();
            grid_step = (*iter)->get("fmatch.step").as<double>();
            i->res_output_coeff = (*iter)->get("fmatch.res_output_coeff").as<int>();
            i->n = i->Spline.GenerateGrid(grid_min, grid_max, grid_step) - 1;
            cout << "Number of splines for the interaction " << i->splineName <<":"<< i->n << endl;
            i->bonded = true;
            i->matr_pos = colms_init;
            i->splineIndex = interaction_number++;
            
            i->result.resize(i->res_output_coeff * (i->n + 1), false);
            i->result.clear();        
            i->error.resize(i->res_output_coeff * (i->n + 1), false);
            i->error.clear();
            i->resSum.resize(i->res_output_coeff * (i->n + 1), false);
            i->resSum.clear();
            i->resSum2.resize(i->res_output_coeff * (i->n + 1), false);
            i->resSum2.clear();
            i->block_res.resize(2*(i->n + 1), false);
            i->del_x_out = (i->Spline.getGridPoint(i->n) - i->Spline.getGridPoint(0)) /
                                    (i->res_output_coeff * (i->n + 1));
            
            //adjust initial matrix dimensions:
            lines_init += i->n + 1;
            colms_init += 2 * (i->n + 1);

            //Add SplineInfo to SplineContainer:
            Splines.push_back( i );
            
            // update bonded interaction counter:
            numBondInt++;
        }
        
        
        // initializing non-bonded interactions
        for (list<Property*>::iterator iter = _nonbonded.begin();
               iter != _nonbonded.end(); ++iter) {
            
            SplineInfo *i = new SplineInfo;
            i->splineName = (*iter)->get("name").value();
            i->type1 = (*iter)->get("type1").value();  // added recently
            i->type2 = (*iter)->get("type2").value();  // !!!
            grid_min = (*iter)->get("fmatch.min").as<double>();
            grid_max = (*iter)->get("fmatch.max").as<double>();
            grid_step = (*iter)->get("fmatch.step").as<double>();
            i->res_output_coeff = (*iter)->get("fmatch.res_output_coeff").as<int>();
            i->n = i->Spline.GenerateGrid(grid_min, grid_max, grid_step) - 1;
            cout << "Number of splines for the interaction " << i->splineName <<":"<< i->n << endl;            
            i->bonded = false;
            i->matr_pos = colms_init;
            i->splineIndex = interaction_number++;
            
            i->result.resize(i->res_output_coeff * (i->n + 1), false);
            i->result.clear();        
            i->error.resize(i->res_output_coeff * (i->n + 1), false);
            i->error.clear();
            i->resSum.resize(i->res_output_coeff * (i->n + 1), false);
            i->resSum.clear();
            i->resSum2.resize(i->res_output_coeff * (i->n + 1), false);
            i->resSum2.clear();
            i->block_res.resize(2*(i->n + 1), false);
            i->del_x_out = (i->Spline.getGridPoint(i->n) - i->Spline.getGridPoint(0)) /
                                    (i->res_output_coeff * (i->n + 1));
            
            //adjust initial matrix dimensions:
            lines_init += i->n + 1;
            colms_init += 2 * (i->n + 1);

            //Add SplineInfo to SplineContainer:
            Splines.push_back( i );
            
        }       
        
/*        // SplineInfo for the first type of bond:
        Bond1.n = Bond1.Spline.GenerateGrid(  0.256, 0.337, 0.005 ) - 1;
        Bond1.bonded = true;
        Bond1.splineIndex = interaction_number++;
        Bond1.splineName = "bond1";
        Bond1.matr_pos = colms_init;
        
        Bond1.res_output_coeff = specify
        Bond1.result.resize(Bond1.res_output_coeff * (Bond1.n + 1), false);
        Bond1.result.clear();
        Bond1.error.resize(Bond1.res_output_coeff * (Bond1.n + 1), false);
        Bond1.error.clear();
        Bond1.resSum.resize(Bond1.res_output_coeff * (Bond1.n + 1), false);
        Bond1.resSum.clear();
        Bond1.resSum2.resize(Bond1.res_output_coeff * (Bond1.n + 1), false);
        Bond1.resSum2.clear();
        Bond2.block_res.resize(2*(Bond2.n + 1), false);
        
        //adjust initial matrix dimensions:
        lines_init += Bond1.n + 1;
        colms_init += 2 * (Bond1.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond1 );
        // update bonded interaction counter:
        numBondInt++;
*/        
        
        // SplineInfo for the second type of bond:
/*        Bond2.n = Bond2.Spline.GenerateGrid( 0.345, 0.395, 0.01) - 1;
        Bond2.bonded = true;
        Bond2.splineIndex = interaction_number++;
        Bond2.splineName = "bond2";
        Bond2.matr_pos = colms_init;
        
        Bond2.res_output_coeff = specify
        Bond2.result.resize(Bond2.res_output_coeff * (Bond2.n + 1), false);
        Bond2.result.clear();
        Bond2.error.resize(Bond2.res_output_coeff * (Bond2.n + 1), false);
        Bond2.error.clear();
        Bond2.resSum.resize(Bond2.res_output_coeff * (Bond2.n + 1), false);
        Bond2.resSum.clear();
        Bond2.resSum2.resize(Bond2.res_output_coeff * (Bond2.n + 1), false);
        Bond2.resSum2.clear();
        Bond2.block_res.resize(2*(Bond2.n + 1), false);
        
        //adjust initial matrix dimensions:
        lines_init += Bond2.n + 1;
        colms_init += 2 * (Bond2.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond2 );
        // update bonded interaction counter:
        numBondInt++;         
*/        
        // SplineInfo for the angle:
/*        Angle1.n = Angle1.Spline.GenerateGrid( 69.8 * 0.0175, 170 * 0.0175, 5 * 0.0175 ) - 1;
        Angle1.bonded = true;
        Angle1.splineIndex = interaction_number++;
        Angle1.splineName = "angle1";
        Angle1.matr_pos = colms_init;

        Angle1.res_output_coeff = specify
        Angle1.result.resize(Angle1.res_output_coeff * (Angle1.n + 1), false);
        Angle1.result.clear();
        Angle1.error.resize(Angle1.res_output_coeff * (Angle1.n + 1), false);
        Angle1.error.clear();
        Angle1.resSum.resize(Angle1.res_output_coeff * (Angle1.n + 1), false);
        Angle1.resSum.clear();
        Angle1.resSum2.resize(Angle1.res_output_coeff * (Angle1.n + 1), false);
        Angle1.resSum2.clear();
        Angle1.block_res.resize(2*(Angle1.n + 1), false);

        //adjust initial matrix dimensions:
        lines_init += Angle1.n + 1;
        colms_init += 2 * (Angle1.n + 1);

        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Angle1 );
        // update bonded interaction counter:
        numBondInt++;        
*/        
/*//===== Non-bonded params ===========================
  //      NB1.n = NB1.Spline.GenerateGrid( 0.279, 1.0, 0.05 ) - 1; // LJ
  //      NB1.n = NB1.Spline.GenerateGrid( 0.229, 0.79, 0.01 ) - 1; //water
        NB1.n = NB1.Spline.GenerateGrid( 0.27, 1.2, 0.02 ) - 1; //methanol
        NB1.bonded = false;
        NB1.splineIndex = interaction_number++;
        NB1.splineName = "non-bonded_1";
        NB1.matr_pos = colms_init;
        
        NB1.res_output_coeff = 5;
        NB1.result.resize(NB1.res_output_coeff * (NB1.n + 1), false);
        NB1.result.clear();
        NB1.error.resize(NB1.res_output_coeff * (NB1.n + 1), false);
        NB1.error.clear();
        NB1.resSum.resize(NB1.res_output_coeff * (NB1.n + 1), false);
        NB1.resSum.clear();
        NB1.resSum2.resize(NB1.res_output_coeff * (NB1.n + 1), false);
        NB1.resSum2.clear();
        NB1.block_res.resize(2*(NB1.n + 1), false);
        
        NB1.del_x_out = (NB1.Spline.getGridPoint(NB1.n) - NB1.Spline.getGridPoint(0)) /
                                    (NB1.res_output_coeff * (NB1.n + 1));
        
        //adjust initial matrix dimensions:
        lines_init += NB1.n + 1;
        colms_init += 2 * (NB1.n + 1);
        
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &NB1 );
 */           
//===================================================*/
  
        N = top->BeadCount(); // Number of beads in topology
        L = 0;                // Initial frame in trajectory  
        excList.CreateExclusions(top);  //exclusion list for non-bonded interactions
        cout << "hey, somebody wants to forcematch!\n";
        
        if (ConstrLeastSQ) { // Constrained Least Squares
            
            // offset, used in EvalConf
            LeastSQOffset = 0;
            
            // B_constr matrix contains continuity conditions for the spline first
            // derivatives.
            B_constr.resize(lines_init, colms_init, false);
            B_constr.clear();
        
            SplineContainer::iterator is;
                
            for(is=Splines.begin(); is != Splines.end(); ++is) {
        
                sfnum = (*is)->n;            
                (*is)->Spline.AddBCToFitMatrix(B_constr, line_cntr, col_cntr);
            
                // update counters
                line_cntr += sfnum + 1;
                col_cntr += 2 * (sfnum + 1);
                
            }
        
            _A.resize( 3*N*N_frames, col_cntr, false); // resize matrix _A
            _b.resize( 3*N*N_frames, false);          // resize vector _b   
            _A.clear();
            _b.clear();        
        
        }
        else {  // Simple Least Squares
            
            // offset, used in EvalConf
            LeastSQOffset = lines_init;            
            
            _A.resize( lines_init + 3*N*N_frames, colms_init, false); // resize matrix _A
            _b.resize( lines_init + 3*N*N_frames, false);          // resize vector _b   
            _A.clear();
            _b.clear();  
            
            SplineContainer::iterator is;
                
            for(is=Splines.begin(); is != Splines.end(); ++is) {
        
                sfnum = (*is)->n;     
                (*is)->Spline.AddBCToFitMatrix(_A, line_cntr, col_cntr);
            
                // update counters
                line_cntr += sfnum + 1;
                col_cntr += 2 * (sfnum + 1);
                
            }            
        }
        
        _x.resize(col_cntr);
        _x.clear(); 
               
    }
    
void CGForceMatching::EndCG() { 
        string force_raw = "_force_raw.dat";
        char file_name[20];
        double accuracy; // accuracy for output. Should be different for bonds and angles.
        
        ofstream out_file;
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
            int &mp = (*is)->matr_pos;
            int &nsf = (*is)->n;
            
            file_name[0] = '\0';
            strcpy(file_name, ((*is)->splineName).c_str() );
            strcat(file_name, force_raw.c_str());
            out_file.open(file_name);
            
            out_file << "# interaction No. " << (*is)->splineIndex << endl;
    
            for (int i = 0; i < (*is)->res_output_coeff * (nsf + 1); i++ ) {
                (*is)->result[i] = (*is)->resSum[i] / BlockNum;
                     if ( i == 23) cout << (*is)->result[i] << endl;
                (*is)->error[i] = sqrt( (*is)->resSum2[i] / BlockNum - (*is)->result[i] * (*is)->result[i]); 
            }

            (*is)->Spline.setSplineData( (*is)->result );
            
            if ( ((*is)->splineName)[0] == 'a' ) accuracy = 0.05;
            else if ( ((*is)->splineName)[0] == 'b' ) accuracy = 0.001;
            else if ( ((*is)->splineName)[0] == 'n' ) accuracy = 0.01; // Quatsch!
            
//            (*is)->Spline.Print(out_file, accuracy);

            // Shitty implementation, think of adding functionality to CubicSpline
            double out_x = (*is)->Spline.getGridPoint(0);
            
            for (int i = 0; i < (*is)->res_output_coeff * (nsf + 1); i++) {
                out_file << out_x << " " <<
                        (-1.0) * (*is)->result[i] << " " << (*is)->error[i] << endl;
                out_x += (*is)->del_x_out;
            }
            
            out_file.close();
        }          
        
        
    }
    

    
void CGForceMatching::EvalConfiguration(Topology *conf, Topology *conf_atom) {
    if (NewNeighbour) { // New neighbour list
        
        SplineContainer::iterator spiter;
        
        for (spiter=Splines.begin(); spiter != Splines.end(); ++spiter) 
        {
            SplineInfo *sinfo = *spiter;
            if (sinfo->bonded) { // bonded interaction
                std::list<Interaction *> interList;
                std::list<Interaction *>::iterator interListIter;
                
                interList = conf->InteractionsInGroup(sinfo->splineName);
                
                for (interListIter=interList.begin(); interListIter!=interList.end();++interListIter) {

                   int beads_in_int = (*interListIter)->BeadCount(); // 2 for bonds, 3 for angles, 4 for dihedrals
                   
                   CubicSpline &SP = sinfo->Spline;

                   int  &mpos = sinfo->matr_pos;
                   int  &nsp = sinfo->n;
               
                   double var = (*interListIter)->EvaluateVar(*conf); // value of bond, angle, or dihedral
                   int i = SP.getInterval(var);   // corresponding spline interval

                   for (int loop = 0; loop < beads_in_int; loop ++) {
                       int ii = (*interListIter)->getBeadId(loop);
                       vec gradient = (*interListIter)->Grad(*conf, loop);
                  
                       SP.AddToFitMatrix(_A, var, 
                               LeastSQOffset + 3*N*L + ii, mpos, gradient.x());
                       SP.AddToFitMatrix(_A, var, 
                               LeastSQOffset + 3*N*L + N + ii, mpos, gradient.y());
                       SP.AddToFitMatrix(_A, var,
                               LeastSQOffset + 3*N*L + 2*N + ii, mpos, gradient.z());                       
                   }
                
                
                }
            }
            else { // non-bonded interaction
                // generate the neighbour list
                NBList NBL;
                NBL.setCutoff(cutoff); // implement different cutoffs for different interactions!
                
                // generate the bead lists
                BeadList beads1, beads2;
                beads1.Generate(*conf, sinfo->type1);
                beads2.Generate(*conf, sinfo->type2); 
                
                // is it same types or different types?
                if(sinfo->type1 == sinfo->type2)
                    NBL.Generate(beads1, &excList);
                else
                    NBL.Generate(beads1, beads2, &excList);
                
                
                NBList::iterator pair_iter;
                // iterate over all pairs
                for(pair_iter = NBL.begin(); pair_iter!=NBL.end();++pair_iter) {
                    int iatom = (*pair_iter)->first->getId();
                    int jatom = (*pair_iter)->second->getId();
                    double var = (*pair_iter)->dist();             
                    vec gradient = (*pair_iter)->r();
                    gradient.normalize();
                    
                    CubicSpline &SP = sinfo->Spline;

                    int  &mpos = sinfo->matr_pos;
                    int  &nsp = sinfo->n;     
                    int i = SP.getInterval(var);
                    
                    // add iatom
                    SP.AddToFitMatrix(_A, var, 
                         LeastSQOffset + 3*N*L + iatom, mpos, gradient.x());
                    SP.AddToFitMatrix(_A, var, 
                         LeastSQOffset + 3*N*L + N + iatom, mpos, gradient.y());
                    SP.AddToFitMatrix(_A, var,
                         LeastSQOffset + 3*N*L + 2*N + iatom, mpos, gradient.z());  
                            
                    // add jatom 
                    SP.AddToFitMatrix(_A, var, 
                         LeastSQOffset + 3*N*L + jatom, mpos, -gradient.x());
                    SP.AddToFitMatrix(_A, var, 
                         LeastSQOffset + 3*N*L + N + jatom, mpos, -gradient.y());
                    SP.AddToFitMatrix(_A, var,
                         LeastSQOffset + 3*N*L + 2*N + jatom, mpos, -gradient.z());                            
                     
                }                
            }
        }
    }       
    else { //Old neighbour list
        
        InteractionContainer &ic = conf->BondedInteractions();
        InteractionContainer::iterator ia;

        // loop for the matrix: (Bonded Interactions)
        for(ia=ic.begin(); ia != ic.end(); ++ia) {
            
               int beads_in_int = (*ia)->BeadCount(); // 2 for bonds, 3 for angles, 4 for dihedrals
                
               int index = (*ia)->getGroupId(); // unique for every interaction type

               CubicSpline &SP = Splines[ index ]->Spline;
               int  &mpos = Splines[ index ]->matr_pos;
               int  &nsp = Splines[ index ]->n;
               
               double var = (*ia)->EvaluateVar(*conf); // value of bond, angle, or dihedral
               int i = SP.getInterval(var);   // corresponding spline interval

               for (int loop = 0; loop < beads_in_int; loop ++) {
                   int ii = (*ia)->getBeadId(loop);
                   vec gradient = (*ia)->Grad(*conf, loop);
                  
                   SP.AddToFitMatrix(_A, var, 
                           LeastSQOffset + 3*N*L + ii, mpos, gradient.x());
                   SP.AddToFitMatrix(_A, var, 
                           LeastSQOffset + 3*N*L + N + ii, mpos, gradient.y());
                   SP.AddToFitMatrix(_A, var,
                           LeastSQOffset + 3*N*L + 2*N + ii, mpos, gradient.z());                
               }
        }
        
        // loop for the matrix: (Nonbonded interactions)
        NeighbourList::container::iterator iter;
        list<int>::iterator excl_iter;
        bool noExcl;
        
        NeighbourList nbl;
        nbl.setCutoff(cutoff);
        nbl.Generate(*conf);
        
        for (int iatom = 0; iatom < N; iatom++) {
            
            if ( excList.GetExclusions(iatom) == NULL ) noExcl = true; 
            else noExcl = false; 
            
            for(iter=nbl.NbList()[iatom]->_neighbours.begin(); iter!=nbl.NbList()[iatom]->_neighbours.end(); iter++){
                int jatom = (*iter)._bead;
                if ( jatom > iatom ) {
                    double var = (*iter)._dist;
                   // cout << var << endl;
                    vec gradient = (*iter)._r;
                    gradient.normalize();
                    
                    
                    if ( !noExcl ) { // iatom has exclusions -> we have to check them
                        list<int> excl_iat = excList.GetExclusions(iatom)->_exclude;
                        for (excl_iter = excl_iat.begin(); excl_iter != excl_iat.end(); ++excl_iter )
                            if ( (*excl_iter) == jatom ) break;
                        if ( excl_iter == excl_iat.end() ) {
                        // iatom and jatom have to be added to matrix

                            int itype = conf->getBead(iatom)->getType()->getId();
                            int jtype = conf->getBead(jatom)->getType()->getId();
                            int int_index = beadType2intType(itype, jtype) + numBondInt;
                            
                            CubicSpline &SP = Splines[ int_index ]->Spline;
                            int  &mpos = Splines[ int_index ]->matr_pos;
                            int  &nsp = Splines[ int_index ]->n;
                            int i = SP.getInterval(var);
                            
                            // add iatom
                            SP.AddToFitMatrix(_A, var, 
                                LeastSQOffset + 3*N*L + iatom, mpos, gradient.x());
                            SP.AddToFitMatrix(_A, var, 
                                LeastSQOffset + 3*N*L + N + iatom, mpos, gradient.y());
                            SP.AddToFitMatrix(_A, var,
                                LeastSQOffset + 3*N*L + 2*N + iatom, mpos, gradient.z());  
                            
                            // add jatom 
                            SP.AddToFitMatrix(_A, var, 
                                LeastSQOffset + 3*N*L + jatom, mpos, -gradient.x());
                            SP.AddToFitMatrix(_A, var, 
                                LeastSQOffset + 3*N*L + N + jatom, mpos, -gradient.y());
                            SP.AddToFitMatrix(_A, var,
                                LeastSQOffset + 3*N*L + 2*N + jatom, mpos, -gradient.z());                            
                   
                        }
                    }
                    else { // iatom has no exclusions. Every neighbor has to be added!
                        // iatom and jatom have to be added to matrix

                        int itype = conf->getBead(iatom)->getType()->getId();
                        int jtype = conf->getBead(jatom)->getType()->getId();
                        int int_index = beadType2intType(itype, jtype) + numBondInt;
                            
                        CubicSpline &SP = Splines[ int_index ]->Spline;
                        int  &mpos = Splines[ int_index ]->matr_pos;
                        int  &nsp = Splines[ int_index ]->n;
                        int i = SP.getInterval(var);
                            
                        // add iatom
                        SP.AddToFitMatrix(_A, var, 
                            LeastSQOffset + 3*N*L + iatom, mpos, gradient.x());
                        SP.AddToFitMatrix(_A, var, 
                            LeastSQOffset + 3*N*L + N + iatom, mpos, gradient.y());
                        SP.AddToFitMatrix(_A, var,
                            LeastSQOffset + 3*N*L + 2*N + iatom, mpos, gradient.z());  
                            
                        // add jatom 
                        SP.AddToFitMatrix(_A, var, 
                            LeastSQOffset + 3*N*L + jatom, mpos, -gradient.x());
                        SP.AddToFitMatrix(_A, var, 
                            LeastSQOffset + 3*N*L + N + jatom, mpos, -gradient.y());
                        SP.AddToFitMatrix(_A, var,
                            LeastSQOffset + 3*N*L + 2*N + jatom, mpos, -gradient.z());                            
                        
                    }
                    
                }
            }
     
        }
    }

        // loop for the forces vector: 
        // hack, chage the Has functions..
        if ( conf->getBead(0)->HasF() ) {
            vec Force(0., 0., 0.);
            for (int iatom = 0; iatom < N; ++iatom) {
                     Force = conf->getBead(iatom)->getF();
                    _b( LeastSQOffset + 3*N*L + iatom) = Force.x();
                    _b( LeastSQOffset + 3*N*L + N+iatom) = Force.y();
                    _b( LeastSQOffset + 3*N*L + 2*N+iatom) = Force.z();
                  //  cout << Force.x() << endl;
            }
        }
        else {
            cout << "ERROR: No forces in configuration!\n" ;   
        }
        L+=1; // update the frame counter
        
        if ( L % N_frames == 0 ) {
            BlockNum++;
            FmatchAccumulateData();
            cout << "Block No" << BlockNum << " done!" << endl;
            L = 0;
            if ( ConstrLeastSQ ) { //Constrained Least Squares
                _A.clear();
                _b.clear();            
            }
            else { // Simple Least Squares
                FmatchAssignMatrixAgain();            
/*                for (int i = line_cntr; i < line_cntr + 3*N*N_frames; i++) {
                    for (int j = 0; j < col_cntr; j++) {
                        _A(i,j)=0.0;
                    }
                }
                _b.clear(); */
            }
        }
    }
    
int CGForceMatching::beadType2intType( int beadType1, int beadType2 ) {
// This function returns the interaction type, knowing the bead types involved.
// The correspondence is established as follows: (case of 4 different bead types)
    
// | interaction | corresponding beads |
// |_____________|_____________________|
// |      0      |       0 - 0         |
// |      1      |       1 - 1         |
// |      2      |       2 - 2         |
// |      3      |       3 - 3         |
// |      4      |       4 - 4         |
// |      5      |       0 - 1         |
// |      6      |       0 - 2         |
// |      7      |       0 - 3         |
// |      8      |       0 - 4         |
// |      9      |       1 - 2          |
// |     10      |       1 - 3         |
// |     11      |       1 - 4         |
// |     12      |       2 - 3         |
// |     13      |       2 - 4         |
// |     14      |       3 - 4         |
// |_____________|_____________________|
    
    int temp, result = 0;

    if ( beadType1 == beadType2 ) return beadType1;
    if ( beadType1 > beadType2 ) {
        temp = beadType1;
        beadType1 = beadType2;
        beadType2 = temp;
    }
    // Now beadType1 < beadType2
    result += beadTypes - 1;
    for ( int i = beadType1 - 1; i >=0; i-- )
        result += beadTypes - 1 - i;
    result += beadType2 - beadType1;
    return result;    
}

    void CGForceMatching::FmatchAccumulateData() {
/*        string force_raw = "_force_raw.dat";
        char file_name[20];
        double accuracy; // accuracy for output. Should be different for bonds and angles.
        
        ofstream out_file;
*/
    //  _x.resize(col_cntr);
        _x.clear(); 
            
        if (ConstrLeastSQ) { // Constrained Least Squares
            
            // Solving linear equations system

            ub::matrix<double> Q; 
            Q.resize(col_cntr, col_cntr, false );
            Q.clear();
        
            ub::matrix<double> A2;  
            A2.resize(_A.size1(), col_cntr / 2, false );
            A2.clear();
                
            ub::matrix<double> Q_k;
            Q_k.resize(col_cntr, col_cntr, false);
            Q_k.clear();
       
            ub::identity_matrix<double> I (col_cntr); 

            ub::vector<double> v; 
            v.resize(col_cntr, false);
            v.clear();        
        
            // To proceed we need to factorize B^T = Q*R. We need matrix Q for further
            // calculations
            B_constr = trans(B_constr);
        
            double* pointer_Bcnstr = & B_constr(0,0); 
        
            gsl_matrix_view B_t 
                = gsl_matrix_view_array (pointer_Bcnstr, col_cntr, line_cntr);
     
            gsl_vector *tau = gsl_vector_alloc (line_cntr); 
    
            gsl_linalg_QR_decomp (&B_t.matrix, tau);   

            // Extraction of Q matrix from tau and B_t, where it is stored in a tricky way.
            Q = I;
              
            for (int k = line_cntr; k > 0 ; k--) {
           
                for (int icout = 0; icout < k - 1; icout++) {
                    v(icout) = 0;
                }
                v(k - 1) = 1.0;

                for (int icout = k; icout < col_cntr; icout++) {
                    v(icout) = gsl_matrix_get(&B_t.matrix, icout, k - 1 );
                }
                double tmp = gsl_vector_get(tau, k - 1 );
                Q_k = I - tmp * outer_prod ( v, v );
                Q = prec_prod(Q, Q_k);
           
            }
    
            
            Q = trans(Q);     
    
            // Calculate _A * Q and store the result in _A
            _A = prec_prod(_A, Q);
    
            // _A = [A1 A2], so A2 is just a block of _A
            for (int iraw = 0; iraw < _A.size1(); iraw++) {
                for (int icol = _A.size2() / 2; icol < _A.size2(); icol++) {
                    A2(iraw, icol - _A.size2() / 2) = _A(iraw, icol);
                }
            }    
        
        
  
            double* pointer_m = & A2(0,0);
            double* pointer_b = & _b(0);        
        
            gsl_matrix_view m
                = gsl_matrix_view_array (pointer_m, A2.size1(), A2.size2() );
    
            gsl_vector_view b
                = gsl_vector_view_array (pointer_b, A2.size1());
    
            gsl_vector *x = gsl_vector_alloc ( A2.size2() );
            gsl_vector *tau2 = gsl_vector_alloc ( A2.size2() );       
            gsl_vector *residual = gsl_vector_alloc ( A2.size1() );
    
            gsl_linalg_QR_decomp (&m.matrix, tau2);
        
            gsl_linalg_QR_lssolve (&m.matrix, tau2, &b.vector, x, residual);
        
            for (int i = 0; i < col_cntr / 2; i++ ) {
                _x[i] = 0.0;
            }  
    
            for (int i = col_cntr / 2; i < col_cntr; i++ ) {
                _x[i] = gsl_vector_get(x, i - col_cntr / 2 );
            }      
    
            // To get the final answer this vector should be multiplied by matrix Q
            _x = prec_prod( Q, _x );       
            
            gsl_vector_free (x);
            gsl_vector_free (tau);
            gsl_vector_free (residual);            
        
        }
        else {  // Simple Least Squares

            double* pointer_m = & _A(0,0);
            double* pointer_b = & _b(0);            

            gsl_matrix_view m
                = gsl_matrix_view_array (pointer_m, _b.size(), col_cntr);
    
            gsl_vector_view b
                = gsl_vector_view_array (pointer_b, _b.size());
    
            gsl_vector *x = gsl_vector_alloc (col_cntr);
            gsl_vector *tau = gsl_vector_alloc (col_cntr);       
            gsl_vector *residual = gsl_vector_alloc (_b.size());
    
            gsl_linalg_QR_decomp (&m.matrix, tau); 
            gsl_linalg_QR_lssolve (&m.matrix, tau, &b.vector, x, residual); 
            
            for (int i =0 ; i < col_cntr; i++) {
                _x(i) = gsl_vector_get(x, i);
            }
            
            gsl_vector_free (x);
            gsl_vector_free (tau);
            gsl_vector_free (residual);                    
            
        }
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
            int &mp = (*is)->matr_pos;
            int &nsf = (*is)->n;
            
            for (int i = 0; i < 2*(nsf + 1); i++ ) {
                (*is)->block_res[i] = _x[ i + mp ];
//                (*is)->resSum[i] += _x[ i + mp ];
//                (*is)->resSum2[i] += _x[ i + mp ] * _x[ i + mp ];
            }
            (*is)->Spline.setSplineData( (*is)->block_res );
            
            double out_x = (*is)->Spline.getGridPoint(0);

            for (int i = 0; i < (*is)->res_output_coeff * (nsf + 1); i++ ) {
                (*is)->resSum[i] += (*is)->Spline.Calculate(out_x);
                    if ( i == 23) cout << (*is)->Spline.Calculate(out_x) << " " << endl;
                (*is)->resSum2[i] += (*is)->Spline.Calculate(out_x) * (*is)->Spline.Calculate(out_x);
                out_x += (*is)->del_x_out;
            }
            
        }         
        
/*        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
            int &mp = (*is)->matr_pos;
            int &nsf = (*is)->n;
            
            file_name[0] = '\0';
            strcpy(file_name, ((*is)->splineName).c_str() );
            strcat(file_name, force_raw.c_str());
            out_file.open(file_name);
            
            out_file << "# interaction No. " << (*is)->splineIndex << endl;
            
            for (int i = 0; i < 2*(nsf + 1); i++ ) {
//               (*is)->result[i] = gsl_vector_get(x, i + mp);
                (*is)->result[i] = _x[ i + mp ];
            }
            //(*is)->Spline.GetResult( & (*is)->result );
            //(*is)->Spline.PrintOutResult();
            (*is)->Spline.setSplineData( (*is)->result );
            
            if ( ((*is)->splineName)[0] == 'a' ) accuracy = 0.05;
            else if ( ((*is)->splineName)[0] == 'b' ) accuracy = 0.001;
            else if ( ((*is)->splineName)[0] == 'n' ) accuracy = 0.01; // Quatsch!
            
            (*is)->Spline.Print(out_file, accuracy);
            
            out_file.close();
        }        
*/        
    }
    
    void CGForceMatching::FmatchAssignMatrixAgain() {
  
            int line_tmp, col_tmp;
            line_tmp = 0;
            col_tmp = 0;
            
            _A.clear();
            _b.clear();  
            
            SplineContainer::iterator is;
                
            for(is=Splines.begin(); is != Splines.end(); ++is) {
        
                int sfnum = (*is)->n;            
                (*is)->Spline.AddBCToFitMatrix(_A, line_tmp, col_tmp);
            
                // update counters
                line_tmp += sfnum + 1;
                col_tmp += 2 * (sfnum + 1);
                
            }             
    }
    
void CGForceMatching::LoadOptions(const string &file) {
    load_property_from_xml(_options, file);   
    _bonded = _options.Select("cg.bonded");
    _nonbonded = _options.Select("cg.non-bonded");
}    
