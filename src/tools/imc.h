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
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
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
        Property *_p;
        HistogramNew _average;
        HistogramNew _current;
        double _min, _max, _step;
        double _norm;
    };
    
    typedef ub::symmetric_matrix<double> group_matrix;
    typedef ub::matrix_range< group_matrix > pair_matrix;
    
    struct pair_t {
        pair_t(interaction_t *i1, interaction_t *i2,
                const pair_matrix &_corr);
                
        interaction_t *_i1;
        interaction_t *_i2;

       pair_matrix _corr;
    };
    
    /// struct to store collected information for groups (e.g. crosscorrelations)
    struct group_t {
        list<interaction_t*> _interactions;
        group_matrix _corr;
    };
    
    /// map ineteractionm-name to interaction
    map<string, interaction_t *> _interactions;
    /// map group-name to group
    map<string, group_t *> _groups;
    vector< pair_t > _pairs;
            
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
    
    // calculate deviation from target vectors
    void CalcDeltaS();
 };
 
 inline Imc::pair_t::pair_t(Imc::interaction_t *i1, Imc::interaction_t *i2,
         const pair_matrix &corr)
         :
  _i1(i1), _i2(i2), _corr(corr)  {}
       


#endif	/* _IMC_H */

