/* 
 * File:   csg_fmatch.h
 * Author: lukyanov
 *
 * Created on June 10, 2009, 4:56 PM
 */

#ifndef _CSG_FMATCH_H
#define	_CSG_FMATCH_H

#include <tools/property.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <tools/cubicspline.h>
#include "cgengine.h"

namespace ub = boost::numeric::ublas;
using namespace std;

class CGForceMatching
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom);
    void EndCG();
    void EvalConfiguration(Topology *conf, Topology *conf_atom = 0);
    void LoadOptions(const string &file);
    
protected:         
  struct SplineInfo {
      SplineInfo(int index, bool bonded_, int matr_pos_, Property *options);
        int n; //number of splines
        int splineIndex; // interaction index for bonded interactions
        bool bonded;     // true for bonded interactions, false for non-bonded
        CubicSpline Spline;
        int matr_pos;    // position in the _A matrix (first coloumn which is occupied with
                         // this particular spline
        int res_output_coeff; // Num_output_points = Num_spline_points * res_output_coeff
        double del_x_out;     // dx for output. Calculated in the code
        
        pair<int, int> beadTypes; // only for non-bonded interactions
        
        
        ub::vector<double> block_res;  // Result of 1 block calculation
        ub::vector<double> result;     // Average over many blocks
        ub::vector<double> error;
        ub::vector<double> resSum;
        ub::vector<double> resSum2;
        
        string splineName;
        string type1, type2; // for non-bonded: types of beads involved
        
        Property *_options;
  };
 
   Property _options;  
  /// list of bonded interactions
   list<Property *> _bonded;
  /// list of non-bonded interactions
   list<Property *> _nonbonded;
      
  typedef vector<SplineInfo *> SplineContainer; 
  SplineContainer Splines;     // please change name to _splines
  
    ub::matrix<double> _A;
  ub::vector<double> _b; // F_ref
  ub::vector<double> _x; // 
  ub::matrix<double> B_constr;    // please change name to _B_constr
  
  // please use better names beginning with _ here, N and L are not very meaningfull
  int L; // counter for frames    
  int N; //number of cg_beads     
  
  // please change to _constr_least_sq or similar
  bool ConstrLeastSQ; // true:  constrained least squares
                      // false: simple least squares
  // please change to _least_sq_offset or similar
  int LeastSQOffset;  // used in EvalConf to distinguish constrained LS and simple LS
  // please change to _nframes or similar
  int N_frames;       // Number of frames used in one Block
  // please change to _nblocks or similar
  int BlockNum;       // current number of Blocks
  
  // please add _ to variable name
  int line_cntr, col_cntr; // counters for lines and coloumns in B_constr 

  
  void FmatchAccumulateData();
  void FmatchAssignMatrixAgain();
  void EvalBonded(Topology *conf, SplineInfo *sinfo);
  void EvalNonbonded(Topology *conf, SplineInfo *sinfo);    
};

#endif	/* _CSG_FMATCH_H */

