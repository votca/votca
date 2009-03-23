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

class Imc
    : public CGObserver
{
public:

    void LoadOptions(const string &file);
    
    void BeginCG(Topology *top, Topology *top_atom);
    
    void EndCG();
    
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);    
    
protected:
    Property _options;
    list<Property *> _bonded;
    list<Property *> _nonbonded;
    
    struct interaction_t {
        HistogramNew _average;
        HistogramNew _current;
        double _min, _max, _step;
    };
    
    struct group_t {
        list<interaction_t*> _interactions;
        vector< ub::matrix<double> > _correlations;
    };
    
    map<string, interaction_t *> _interactions;
    map<string, group_t *> _groups;
            
    int _nframes;
    
    interaction_t *AddInteraction(Property *p);
    
    void DoNonbonded(Topology *top);
    group_t *getGroup(const string &name);
    
    void InitializeGroups();
    void UpdateCorrelations();

    void Correlate(ub::matrix<double> &corr, ub::vector<double> &v1,
        ub::vector<double> &v2);
 };


#endif	/* _IMC_H */

