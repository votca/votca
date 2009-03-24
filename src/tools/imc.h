/* 
 * File:   imc.h
 * Author: ruehle
 *
 * Created on March 10, 2009, 3:42 PM
 */

#ifndef _IMC_H
#define	_IMC_H

#include <cgengine.h>
#include <tools/property.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <tools/histogramnew.h>

namespace ub = boost::numeric::ublas;
using namespace std;

/**
 * \brief class to calculate distribution functions and cross correlations for inverse monte carlo
 * 
 * This class calculates distribution functions as well as cross-correlations for specific groups of
 * interactions based on a given trajectory.
 * 
 */
class Imc
    : public CGObserver
{
public:

    /// load cg definitions file
    void LoadOptions(const string &file);
    
    /// begin coarse graining a trajectory
    void BeginCG(Topology *top, Topology *top_atom);
    
    /// end coarse graining a trajectory
    void EndCG();
    
    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);    
    
protected:
    /// the options parsed from cg definition file
    Property _options;
    
    /// list of bonded interactions
    list<Property *> _bonded;
    /// list of non-bonded interactions
    list<Property *> _nonbonded;
    
    /// struct to store collected information for interactions
    struct interaction_t {
        HistogramNew _average;
        HistogramNew _current;
        double _min, _max, _step;
    };
    
    /// struct to store collected information for groups (e.g. crosscorrelations)
    struct group_t {
        list<interaction_t*> _interactions;
        vector< ub::matrix<double> > _correlations;
    };
    
    /// map ineteractionm-name to interaction
    map<string, interaction_t *> _interactions;
    /// map group-name to group
    map<string, group_t *> _groups;
            
    int _nframes;
    
    /// create a new interaction entry based on given options
    interaction_t *AddInteraction(Property *p);
    
    /// process non-bonded interactions for given frame
    void DoNonbonded(Topology *top);
    
    /// get group by name, creates one if it doesn't exist
    group_t *getGroup(const string &name);
    
    /// initializes the group structs after interactions were added
    void InitializeGroups();
    /// update the correlations after interations were processed
    void UpdateCorrelations();

    /// update the correlations after interations were processed
    void CalcMatrix();
 };


#endif	/* _IMC_H */

