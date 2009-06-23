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
    // _A*_x = _b
    
  ub::matrix<double> _A;
  ub::vector<double> _b; // F_ref
  ub::vector<double> _x; // 
  ub::matrix<double> B_constr;
    // _A(i, j) = 10;
    // _A.resize(n, m, true);
  int beadTypes; // number of cg bead types in the system
  int numBondInt; // number of bonded interaction types
  int L; // counter for frames
  int N; //number of cg_beads
  bool ConstrLeastSQ; // true:  constrained least squares
                      // false: simple least squares
  int LeastSQOffset;  // used in EvalConf to distinguish constrained LS and simple LS
  int N_frames;       // Number of frames used in one Block
  int BlockNum;       // current number of Blocks
  
  int line_cntr, col_cntr; // counters for lines and coloumns in B_constr 
  
  struct SplineInfo {
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
  };
 // SplineInfo Bond1;
//  SplineInfo Bond2;
//  SplineInfo Angle1;
//  SplineInfo Angle2;
//  SplineInfo Angle3;
  SplineInfo NB1;
  
// NEW STUFF
  /// the options parsed from cg definition file
   Property _options;  
  
  /// list of bonded interactions
   list<Property *> _bonded;
  /// list of non-bonded interactions
   list<Property *> _nonbonded;
   
// END NEW STUFF
 
  ExclusionList excList;  // exclusion list for non-bonded interactions
  
  typedef vector<SplineInfo *> SplineContainer;
  SplineContainer Splines;
  
  int beadType2intType ( int beadType1, int beadType2 );
  void FmatchAccumulateData();
  void FmatchAssignMatrixAgain();
  
};

#endif	/* _CSG_FMATCH_H */

