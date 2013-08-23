
#include <votca/csg/potentialfunctions/potentialfunctioncbspl.h>

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const string& name_,const int nlam_,const string & core_,
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

  if( core_ != "c12" && core_ != "extrapolate" )
    throw std::runtime_error("Repulsive core type "+core_+ "selected for potential "+_name+" "
                             "does not exist! \n "
                             "Choose either c12 or extrapolate. ");

  _core = core_;

  int nknots;
  double rstart;


  if( _core == "c12" ){

    nknots = _lam.size() - 1;
    rstart = _min;
    _nexcl = 0;

  } else {

    nknots = _lam.size();
    rstart = 0.0;

  }

  _nbreak = nknots - 2;

  _dr = ( _cut_off - rstart )/( double (_nbreak - 1) );

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak.resize(nknots,false);
  _rbreak.clear();

  for( int i = 0; i < nknots; i++)
    _rbreak(i) = rstart + i*_dr;

  if( _core == "extrapolate" ) {

    // exclude knots corresponding to r <= _min
    _nexcl = min( int( ( _min )/_dr ), _nbreak - 2 ) + 1;

    // account for finite numerical division of _min/_dr
    // e.g. 0.24/0.02 may result in 11.99999999999999
    if( _rbreak(_nexcl) == _min ) _nexcl++;

  }

  // fixing last 4 knots to zeros is reasonable
  _ncutcoeff = 4;

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
    for( int i = 0; i < _lam.size() - _ncutcoeff; i++)
      _lam(i) = param.y(i);

  }

}

void PotentialFunctionCBSPL::SaveParam(const string& filename){

  if( _core == "extrapolate" )
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
  if( _core == "extrapolate" )
    extrapolExclParam();

  PotentialFunction::SavePotTab(filename,step,rmin,rcut);

}

void PotentialFunctionCBSPL::SavePotTab(const string& filename,
                                        const double step) {

  if( _core == "extrapolate" )
    extrapolExclParam();

  PotentialFunction::SavePotTab(filename,step);

}

void PotentialFunctionCBSPL::extrapolExclParam(){

  double u0 = _lam(_nexcl);
  double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
    (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
  double r0 = _rbreak(_nexcl);

  // preferred extrapolation function is exponential
  // exponential extrapolation is useful only if u0 is positive
  if( u0 > 0.0 ){

    cout << "Using exponential function to extrapolate "
      "spline knots in the repulsive core"<< endl;

    // u(r) = a * exp( b * r)
    // a = u0 * exp ( - m * r0/u0 )
    // b = m/u0
    // m = (u1-u0)/(r1-r0)
    double a = u0 * exp(-m * r0 / u0);
    double b = m / u0;
    for (int i = 0; i < _nexcl; i++)
      _lam(i)  = a * exp(b * _rbreak(i));

  } else if( m < 0.0 ) // try linear extrapolation (physical, i.e. repulsive, only when slope is negative)
    {

      // u(r) = ar + b
      // a = m
      // b = - m*r0 + u0
      // m = (u1-u0)/(r1-r0)

      double a = m;
      double b = -1.0*m*r0 + u0;
      for (int i = 0; i < _nexcl; i++)
        _lam(i) = a*_rbreak(i) + b;

    } else // extrapolation results in unphysical, i.e., non-repulsive, core
    {

      throw std::runtime_error("In potential "+_name+": extrapolation results in unphysical potential core!\n"
                               "Check if rmin value for cbspl is too large.\n");

    }

}

void PotentialFunctionCBSPL::setOptParam(const int i, const double val){

  _lam( i + _nexcl ) = val;

}

double PotentialFunctionCBSPL::getOptParam(const int i) const{

  return _lam( i + _nexcl );

}

double PotentialFunctionCBSPL::CalculateF (const double r) const {

  if( r <= _cut_off){

    double u = 0.0;

    ub::vector<double> R;
    ub::vector<double> B;
    R.resize(4,false); R.clear();
    B.resize(4,false); B.clear();

    int indx;
    double rk;

    if( _core == "c12" ){

      u += _lam(0)/pow(r,12);

      indx = min( int( (r - _min)/_dr ), _nbreak-2);

      if( indx >= 0){ // r >= _min

        indx += 1;
        rk = _min + (indx-1)*_dr;

      }

    }else{

      indx = min( int( r /_dr ) , _nbreak-2 );
      rk = indx*_dr;

    }

    if( indx >= 0 ){

      double t = ( r - rk)/_dr;

      R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;
      ub::vector<double> RM = ub::prod(R,_M);
      B(0) = _lam(indx); B(1) = _lam(indx+1); B(2) = _lam(indx+2);
      B(3) = _lam(indx+3);

      u += ub::inner_prod(B,RM);

    }

    return u;

  } else
    return 0.0;

}

// calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const int i, const double r) const{

  // since first _nexcl parameters are not optimized for stability reasons
  //i = i + _nexcl;

  if ( r <= _cut_off ) {

    int indx, i_indx;
    double rk;

    if( _core == "c12" ){

      if( i == 0 )
        return 1.0/pow(r,12.0);

      indx = min( int( (r - _min)/_dr ), _nbreak-2);

      if( indx >= 0){ // r >= _min

        indx += 1;
        rk = _min + (indx-1)*_dr;

      }

    }else{

      indx = min( int( ( r )/_dr ), _nbreak-2 );
      rk = indx*_dr;

    }

    if( indx >= 0 ){

      i_indx = i + _nexcl;

      if ( i_indx >= indx && i_indx <= indx+3 ){

        ub::vector<double> R;
        R.resize(4,false); R.clear();

        double t = ( r - rk)/_dr;

        R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

        ub::vector<double> RM = ub::prod(R,_M);

        return RM(i_indx-indx);

      }else
        return 0.0;

    } else
      return 0.0;

  } else
    return 0;

}

// calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const int i, const int j,
                                            const double r) const {

  return 0.0;

}
