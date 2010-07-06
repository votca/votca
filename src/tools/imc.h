/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _IMC_H
#define	_IMC_H

#include <cgengine.h>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <votca/tools/histogramnew.h>

namespace votca { namespace csg {
using namespace votca::tools;

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
    Imc();
    ~Imc();
    
    /// load cg definitions file
    void LoadOptions(const string &file);
    
    /// begin coarse graining a trajectory
    void BeginCG(Topology *top, Topology *top_atom);
    
    /// end coarse graining a trajectory
    void EndCG();
    
    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);    
    
    
    void WriteEvery(int write_every) { _write_every = write_every; }
    void DoBlocks(bool do_blocks) { _do_blocks = do_blocks; }
    void DoImc(bool do_imc) { _do_imc = do_imc; }
    
protected:
    
    typedef ub::matrix<double> group_matrix;
    typedef ub::matrix_range< group_matrix > pair_matrix;
    
    /// struct to store collected information for interactions
    struct interaction_t {
        Property *_p;
        HistogramNew _average;
        HistogramNew _current;
        double _min, _max, _step;
        double _norm;
        bool _is_bonded;
    };
    
    // a pair of interactions which are correlated
    struct pair_t {
        pair_t(interaction_t *i1, interaction_t *i2,
                int offset_i, int offset_j, const pair_matrix &corr);                
        interaction_t *_i1;
        interaction_t *_i2;
        pair_matrix _corr;
        int _offset_i, _offset_j;
    };
    
    /// struct to store collected information for groups (e.g. crosscorrelations)
    struct group_t {
        list<interaction_t*> _interactions;
        group_matrix _corr;
        vector<pair_t> _pairs;
    };
    
    
    /// the options parsed from cg definition file
    Property _options;
    // we want to write out every so many frames
    int _write_every;
    // we want do do block averaging -> clear averagings every write out
    bool _do_blocks;
    // calculate the inverse monte carlos parameters (cross correlations)
    bool _do_imc;

    // number of frames we processed
    int _nframes;
    int _nblock;

    /// list of bonded interactions
    list<Property *> _bonded;
    /// list of non-bonded interactions
    list<Property *> _nonbonded;        
    
    /// map ineteractionm-name to interaction
    map<string, interaction_t *> _interactions;
    /// map group-name to group
    map<string, group_t *> _groups;
                    
    /// create a new interaction entry based on given options
    interaction_t *AddInteraction(Property *p);
        
    /// get group by name, creates one if it doesn't exist
    group_t *getGroup(const string &name);
    
    /// initializes the group structs after interactions were added
    void InitializeGroups();    

    /// process non-bonded interactions for given frame
    void DoNonbonded(Topology *top);
    /// process bonded interactions for given frame
    void DoBonded(Topology *top);
    /// update the correlations after interations were processed
    void DoCorrelations();
    
    void WriteDist(const string &suffix="");
    void WriteIMCData(const string &suffix="");
    void WriteIMCBlock(const string &suffix);

    void CalcDeltaS(interaction_t *interaction, 
                    ub::vector_range< ub::vector<double> > &dS);
    
    void ClearAverages();
 };
 
 inline Imc::pair_t::pair_t(Imc::interaction_t *i1, Imc::interaction_t *i2,
         int offset_i, int offset_j, const pair_matrix &corr)
         :
  _i1(i1), _i2(i2), _offset_i(offset_i), _offset_j(offset_j), _corr(corr)  {}

}}

#endif	/* _IMC_H */

