

#include <votca/tools/table.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

using namespace std;
namespace po = boost::program_options;
using namespace votca::tools;

/**
 *  *** Fit a cubic B-spline to given data ***
 *
 */

void help_text(void)
{
  cout << "This program fits a cubic B-spline to input x-y data, and \n"
    "outputs the knot values, break points, and fitted cubic B-spline function.\n\n";
}

void check_option(po::options_description &desc, po::variables_map &vm, const string &option)
{
  if(!vm.count(option)) {
    cout << "votca_fitcbspl \n\n";
    cout << desc << endl << "parameter " << option << " is not specified\n";
    exit(1);
  }
}

int main(int argc, char** argv)
{


  // read in the inputs
  string data_file;
  int ncoeff;

  // Declare the supported options.
  po::options_description desc("Allowed options");
    
  // let cg_engine add some program options
  desc.add_options()
    ("data", po::value<string>(&data_file), "target x-y data file")
    ("nknots", po::value<int>(&ncoeff), "number of knots for cubic B-spline\n" 
     "( nbreak = nknots -2 )")
    ("help", "produce this help message");
    
  // now read in the command line
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  }
  catch(po::error err) {
    cout << "error parsing command line: " << err.what() << endl;
    return -1;
  }

  // does the user want help?
  if (vm.count("help")) {
    help_text();
    cout << desc << endl;
    return 0;
  }
  
  // check if all required options are provided
  check_option(desc, vm, "data");
  check_option(desc, vm, "nknots");

  // read in the input x-y data
  Table data;
  data.Load(data_file);

  int ndata = data.size();

  // determine x limits
  double x_min = data.x(0);
  double x_max = data.x(ndata-1);

  // for cubic B-spline nbreak = ncoeff - 2
  int nbreak =  ncoeff - 2;

  // break interval size
  double dx = (x_max - x_min )/( double (nbreak - 1) );

  // cubic B-spline computation matrix M
  ub::matrix<double> M;
  M.resize(4,4,false);
  M.clear();
  M(0,0) =  1.0; M(0,1) =  4.0; M(0,2) =  1.0; M(0,3) = 0.0;
  M(1,0) = -3.0; M(1,1) =  0.0; M(1,2) =  3.0; M(1,3) = 0.0;
  M(2,0) =  3.0; M(2,1) = -6.0; M(2,2) =  3.0; M(2,3) = 0.0;
  M(3,0) = -1.0; M(3,1) =  3.0; M(3,2) = -3.0; M(3,3) = 1.0;
  M /= 6.0;

  // gsl stuff
  gsl_vector *y,*coeff;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq;

  y = gsl_vector_alloc(ndata);
  coeff = gsl_vector_alloc(ncoeff);
  X = gsl_matrix_alloc(ndata, ncoeff);
  cov = gsl_matrix_alloc(ncoeff, ncoeff);
  mw = gsl_multifit_linear_alloc(ndata, ncoeff);

  gsl_matrix_set_zero(X);

  // construct the fit matrix X
  for (int i = 0; i < ndata; i++){

    double xi = data.x(i);
    double indx = min( int( ( xi - x_min )/dx ), nbreak-2 );
    double xk   = x_min + indx*dx;
    double t = ( xi - xk)/dx;

    ub::vector<double> R;
    R.resize(4,false); R.clear();
    R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

    ub::vector<double> RM = ub::prod(R,M);

    // fill in row i of X
    int k = 0;
    for (int j = indx; j < indx + 4; ++j){

      gsl_matrix_set(X, i, j, RM(k));
      k++;

    }
    // also store pot values in gsl_vector y
    gsl_vector_set(y, i, data.y(i));

  }

  /* do the fit */
  gsl_multifit_linear(X, y, coeff, cov, &chisq, mw);

  // write out knot values in knots.dat
  ofstream fl_out;
  string out_file = "knots.dat";

  try {
    fl_out.open(out_file.c_str());
    if(!fl_out.is_open())
      throw std::runtime_error("can't open " + out_file);

    for( int i = 0; i < ncoeff; i++ ) {
      
      double xi =  x_min + i*dx;
      fl_out<< xi << "\t" <<  gsl_vector_get(coeff,i) << endl;

    }
    fl_out.close();
  }
  catch(std::exception &error) {
    cerr << "An error occured!" << endl << error.what() << endl;
  }

  // write out fitted cubic B-spline function values in cbsplfit.dat
  out_file = "cbsplfit.dat";

  try {
    fl_out.open(out_file.c_str());
    if(!fl_out.is_open())
      throw std::runtime_error("can't open " + out_file);

    for( int i = 0; i < ndata; i++ ) {
      
      double xi =  data.x(i);
      double indx = min( int( ( xi - x_min )/dx ), nbreak-2 );
      double xk   = x_min + indx*dx;
      double t = ( xi - xk)/dx;

      ub::vector<double> R;
      R.resize(4,false); R.clear();
      R(0) = 1.0; R(1) = t; R(2) = t*t; R(3) = t*t*t;

      ub::vector<double> RM = ub::prod(R,M);

      ub::vector<double> B;
      B.resize(4,false); B.clear();
      B(0) = gsl_vector_get(coeff,indx);
      B(1) = gsl_vector_get(coeff,indx+1); 
      B(2) = gsl_vector_get(coeff,indx+2);
      B(3) = gsl_vector_get(coeff,indx+3);

      double u = ub::inner_prod(B,RM);
      
      fl_out<< xi << "\t" <<  u << endl;

    }
    fl_out.close();
  }
  catch(std::exception &error) {
    cerr << "An error occured!" << endl << error.what() << endl;
  }

  cout << "Cubic B-spline fitted to the input data" << endl;
  cout << "Fit quality:" << endl;
  cout << "\t \t" << "chisq = "<< chisq << endl;
  cout << "Output in files \"knots.dat\" and \"cbsplfit.dat\""<<endl;
  

  // free gsl specific pointers
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(coeff);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

}
