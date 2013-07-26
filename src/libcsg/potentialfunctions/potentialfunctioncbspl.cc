
#include <votca/csg/potentialfunctions/potentialfunctioncbspl.h>

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const string& name_,const int nlam_,
                                               const int ncutcoeff_, const string& extrapol_type_,
                                               const double min_, const double max_) :
  PotentialFunction(name_,nlam_,min_,max_) {

  /* Here nlam_ is the total number of coeff values that are to be optimized
   * To ensure that potential and force go to zero smoothly near cut-off,
   * as suggested in Ref. PCCP, 11, 1901, 2009, coeff values leading up to
   * cut-off and beyond take a value of zero.
   *
   * Since region less than rmin is not sampled sufficiently for stability
   * first _nexcl coefficients are not optimized instead their values are
   * extrapolated from first statistically significant knot values near rmin
   */

  /* determine minimum for B-spline from CG-MD rdf */
  Table dist;
  string filename = _name + ".dist.cur";
  dist.Load(filename);
  double new_min;

  for( int i = 0; i < dist.size(); i++){

    if(dist.y(i) > 1.0e-3){
      new_min = dist.x(i);
      break;
    }

  }

  if(new_min > _min)
    _min = new_min;

  _extrapol_type = extrapol_type_;
  // number of break points = _lam.size() - 2
  _nbreak = _lam.size() - 2;

  _dr = ( _cut_off )/( double (_nbreak - 1) );

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak.resize(_lam.size(),false);

  _rbreak.clear();

  for( int i = 0; i < _lam.size(); i++)
    _rbreak(i) = i*_dr;

  // exclude knots corresponding to r <= _min
  _nexcl = min( int( ( _min )/_dr ), _nbreak - 2 ) + 1;

  // account for finite numerical division of _min/_dr
  // e.g. 0.24/0.02 may result in 11.99999999999999
  if( _rbreak(_nexcl) == _min ) _nexcl++;

  _ncutcoeff = ncutcoeff_;

  _M.resize(4,4,false);
  _M.clear();
  _M(0,0) =  1.0; _M(0,1) =  4.0; _M(0,2) =  1.0; _M(0,3) = 0.0;
  _M(1,0) = -3.0; _M(1,1) =  0.0; _M(1,2) =  3.0; _M(1,3) = 0.0;
  _M(2,0) =  3.0; _M(2,1) = -6.0; _M(2,2) =  3.0; _M(2,3) = 0.0;
  _M(3,0) = -1.0; _M(3,1) =  3.0; _M(3,2) = -3.0; _M(3,3) = 1.0;
  _M /= 6.0;

}

int PotentialFunctionCBSPL::getOptParamSize() const {

  return _lam.size() - _nexcl - _ncutcoeff;

}

void PotentialFunctionCBSPL::setParam(string filename) {

  Table param;
  param.Load(filename);
  _lam.clear();

  if( param.size() != _lam.size()) {

    throw std::runtime_error("In potential "+_name+": parameters size mismatch!\n"
                             "Check input parameter file \""
                             + filename + "\" \nThere should be "
                             + boost::lexical_cast<string>( _lam.size() ) + " parameters");
  } else {
    // force last _ncutcoeff to be zero
    for( int i = 0; i < _lam.size() - _ncutcoeff; i++){

      _rbreak(i) = param.x(i);
      _lam(i) = param.y(i);

    }

  }

}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  extrapolExclParam();

  Table param;
  param.SetHasYErr(false);
  param.resize(_lam.size(), false);

  // write extrapolated knots with flag 'o'
  for (int i = 0; i < _nexcl; i++)
    param.set(i, _rbreak(i), _lam(i), 'o');

  for (int i = _nexcl; i < _lam.size(); i++)
    param.set(i, _rbreak(i), _lam(i), 'i');

  param.Save(filename);

}

void PotentialFunctionCBSPL::SavePotTab(const string& filename,
                                        const double step, const double rmin, const double rcut) {

  extrapolExclParam();
  PotentialFunction::SavePotTab(filename,step,rmin,rcut);

}

void PotentialFunctionCBSPL::SavePotTab(const string& filename,
                                        const double step) {

  extrapolExclParam();
  PotentialFunction::SavePotTab(filename,step);

}

void PotentialFunctionCBSPL::extrapolExclParam(){

  if( _extrapol_type == "linear") {
    // extrapolate first _nexcl knot values using linear extrapolation
    // u(r) = ar + b
    // a = m
    // b = - m*r0 + u0
    // m = (u1-u0)/(r1-r0)
    double u0 = _lam(_nexcl);
    double r0 = _rbreak(_nexcl);
    double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
      (_rbreak(_nexcl + 1) - _rbreak(_nexcl));

    // check for positive slope
    if( m >= 0.0 )
      throw std::runtime_error("In potential "+_name+": min r value for cbspl is too large,\n"
                               "for linear extrapolation, choose min r such that potential slope at min < 0,\n"
                               "else extrapolated knot values in "
                               "the repulsive core would be negative or zero.\n");

    double a = m;
    double b = -1.0*m*r0 + u0;
    for (int i = 0; i < _nexcl; i++)
      _lam(i) = a*_rbreak(i) + b;

  } else if (_extrapol_type=="exponential"){
    // extrapolate first _nexcl knot values using exponential extrapolation
    // u(r) = a * exp( b * r)
    // a = u0 * exp ( - m * r0/u0 )
    // b = m/u0
    // m = (u1-u0)/(r1-r0)
    double u0 = _lam(_nexcl);

    if( u0 <= 0.0 )
      throw std::runtime_error("In potential "+_name+": min r value for cbspl is too large,\n"
                               "for exponential extrapolation, choose min r such that knot value at (rmin+dr) > 0,\n"
                               "else extrapolated knot values in "
                               "the repulsive core would be negative or zero.\n");

    double r0 = _rbreak(_nexcl);
    double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
      (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
    double a = u0 * exp(-m * r0 / u0);
    double b = m / u0;
    for (int i = 0; i < _nexcl; i++)
      _lam(i)  = a * exp(b * _rbreak(i));

  } else
    throw std::runtime_error("Extrapolation method  \""
                             + _extrapol_type + "\" selected for \""
                             + _name + "\" is not available yet.\n"
                             + "Please specify either \"linear or exponential\" "
                             + " in options file.");

}

void PotentialFunctionCBSPL::setOptParam(const int i, const double val){

  _lam( i + _nexcl ) = val;

}

double PotentialFunctionCBSPL::getOptParam(const int i) const{

  return _lam( i + _nexcl );

}

double PotentialFunctionCBSPL::CalculateF (const double r) const {

  if( r <= _cut_off){

    ub::vector<double> R;
    ub::vector<double> B;
    R.resize(4,false); R.clear();
    B.resize(4,false); B.clear();

    int indx = min( int( r /_dr ) , _nbreak-2 );
    double rk = indx*_dr;
    double t = ( r - rk)/_dr;

    R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;
    ub::vector<double> RM = ub::prod(R,_M);
    B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
    B(3) = _lam(indx+3);

    return ub::inner_prod(B,RM);

  } else
    return 0.0;

}

// calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const int i, const double r) const{

  // since first _nexcl parameters are not optimized for stability reasons
  //i = i + _nexcl;

  if ( r <= _cut_off ) {

    int indx = min( int( ( r )/_dr ), _nbreak-2 );

    if ( i + _nexcl >= indx && i + _nexcl <= indx+3 ){

      ub::vector<double> R;
      R.resize(4,false); R.clear();

      double rk = indx*_dr;
      double t = ( r - rk)/_dr;

      R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

      ub::vector<double> RM = ub::prod(R,_M);

      return RM(i + _nexcl-indx);

    }else
      return 0.0;

  } else
    return 0;

}

// calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const int i, const int j,
                                            const double r) const {

  // for cubic B-SPlines D2F is zero for all lamdas
  return 0.0;

}
