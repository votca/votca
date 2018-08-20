/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/csg/csgapplication.h>
#include <votca/tools/property.h>
#include <votca/tools/histogramnew.h>
#include <votca/tools/average.h>

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
{
public:
    Imc();
    ~Imc();

    void Initialize(void);

    /// load cg definitions file
    void LoadOptions(const string &file);
    
    /// begin coarse graining a trajectory
    void BeginEvaluate(Topology *top, Topology *top_atom);
    
    /// end coarse graining a trajectory
    void EndEvaluate();
    
    void BlockLength(int length) { _block_length = length; }
    void DoImc(bool do_imc) { _do_imc = do_imc; }
    void Extension(string ext) { _extension = ext; }
    
protected:
    Average<double> _avg_vol;
    
    typedef Eigen::MatrixXd group_matrix;
    typedef Eigen::Block<group_matrix > pair_matrix;
    
    /// struct to store collected information for interactions
    struct interaction_t {
        int _index;
        Property *_p;
        HistogramNew _average;
        HistogramNew _average_force;
        /*HistogramNew _average_force_perp;
        HistogramNew _average_force_perp_dot;        
        HistogramNew _average_force_perp_x;
        HistogramNew _average_force_perp_y;
        HistogramNew _average_force_perp_z;*/
        double _min, _max, _step;
        double _norm;
        double _cut;
        bool _is_bonded;
        bool _threebody;
        bool _force;
    };
    
    // a pair of interactions which are correlated
    struct pair_t {
        interaction_t *_i1;
        interaction_t *_i2;
        int _offset_i, _offset_j;
        pair_matrix _corr;
        pair_t(interaction_t *i1, interaction_t *i2,
                int offset_i, int offset_j, const pair_matrix &corr);                
    };
    
    /// struct to store collected information for groups (e.g. crosscorrelations)
    struct group_t {
        list<interaction_t*> _interactions;
        group_matrix _corr;
        vector<pair_t> _pairs;
    };
    
    
    /// the options parsed from cg definition file
    Property _options;
    // length of the block to write out and averages are clear after every write
    int _block_length;
    // calculate the inverse monte carlos parameters (cross correlations)
    bool _do_imc;

    // file extension for the distributions
    string _extension;

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

    
    void WriteDist(const string &suffix="");
    void WriteIMCData(const string &suffix="");
    void WriteIMCBlock(const string &suffix);

    void CalcDeltaS(interaction_t *interaction, 
                    Eigen::VectorBlock< Eigen::VectorXd > &dS);
    
    void ClearAverages();

    class Worker : public CsgApplication::Worker
    {
    public:

        vector<HistogramNew> _current_hists;
        vector<HistogramNew> _current_hists_force;
        /*vector<HistogramNew> _current_hists_force_perp;
        vector<HistogramNew> _current_hists_force_perp_dot;        
        vector<HistogramNew> _current_hists_force_perp_x;
        vector<HistogramNew> _current_hists_force_perp_y;
        vector<HistogramNew> _current_hists_force_perp_z;*/
        Imc *_imc;
        double _cur_vol;

        /// evaluate current conformation
        void EvalConfiguration(Topology *top, Topology *top_atom);
        /// process non-bonded interactions for given frame
        void DoNonbonded(Topology *top);
        /// process bonded interactions for given frame
        void DoBonded(Topology *top);
    };
/// update the correlations after interations were processed
    void DoCorrelations(Imc::Worker *worker);

    bool _processed_some_frames;
    
public:
    CsgApplication::Worker *ForkWorker();
    void MergeWorker(CsgApplication::Worker *worker);
 };
 
 inline Imc::pair_t::pair_t(Imc::interaction_t *i1, Imc::interaction_t *i2,
         int offset_i, int offset_j, const pair_matrix &corr)
         :
  _i1(i1), _i2(i2), _offset_i(offset_i), _offset_j(offset_j), _corr(corr)  {}

}}

#endif	/* _IMC_H */

