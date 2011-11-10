/*
 * File:   doublesum_f.h
 * Author: mashaya1
 *
 * Created on November 8, 2011, 11:48 PM
 */

/* This contains function which computes the double sum of
 * any pair function between beads 
 * i.e. \sum_i \sum_j f(r_ij)
 */

#ifndef DOUBLESUM_F_H
#define	DOUBLESUM_F_H

#include <votca/csg/beadlist.h>
#include "potfunctions.h"
#include "typedefpotfun.h"
#include <votca/csg/nblistgrid.h>

using namespace std;

// potential function class member function pointer
typedef  double (PotFunction::*PotMemFn)(double );

// Macro to call potential function class member function through its pointer
#define CALL_MEMBER_FN(object,ptrToMember)  ((object)->*(ptrToMember))

// from beadlist for two different bead types
double DoubleSum_F(BeadList beads1, BeadList beads2,
                        PotFunction* u, PotMemFn f) {

    // caller function should make sure that beads1 and beads2 has non-zero size

    double sumF = 0.0;

    // generate neighbor list
    NBList *nb;
    nb->setCutoff(u->getCutOff());
    nb->Generate(beads1, beads2, true);

    // compute double sum by iterating over all pairs
    NBList::iterator pair_iter;
    for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
        // assuming that nb->dist() takes care of pbc
        double dist = (*pair_iter)->dist();
        sumF += CALL_MEMBER_FN(u,f)(dist);
    }

    delete nb;
    
    return sumF;

}

// from beadlist for the same bead types
double DoubleSum_F(BeadList beads, PotFunction* u, PotMemFn f) {

    // caller function should make sure that beads has non-zero size

    double sumF = 0.0;

    // generate neighbor list
    NBList *nb;
    nb->setCutoff(u->getCutOff());
    nb->Generate(beads, true);

    // compute double sum \sum_i \sum_j f(r_ij)
    // by iterating over all pairs
    NBList::iterator pair_iter;
    for (pair_iter = nb->begin(); pair_iter != nb->end(); ++pair_iter) {
        // assuming that nb->dist() takes care of pbc
        double dist = (*pair_iter)->dist();
        sumF += CALL_MEMBER_FN(u,f)(dist);
    }

    delete nb;

    return sumF;

}

// using histogram
double DoubleSum_F(Table &hist,int nbeads, bool samebeads,
        PotFunction* u, PotMemFn f){

    /*
     * with histogram data
     * if histogram is for pair b1-b2 then nbeads = nbeads of b1 and vice versa
     * 
     * \sum_i \sum_j f(r_ij) = nbeads*\sum_k hist(r_k)f(r_k)
     *
     * if histogram is for pair b1-b1 then nbeads = nbeads of b1
     *
     * \sum_i \sum_j f(r_ij) = (nbeads*\sum_k hist(r_k)f(r_k))/2.0 to avoid
     *                                                              duplication
     */

    double sumF = 0.0;

    // compute sum
    for(int i = 0; i < hist.size(); i++){
        sumF += hist.y(i)*CALL_MEMBER_FN(u,f)(hist.x(i));
    }

    if( samebeads ){
        return nbeads*sumF/2.0;
    }else {
        return nbeads*sumF;
    }
   
}

#endif	/* DOUBLESUM_F_H */