// 
// File:   csg_nemat.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:35 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <tools/cubicspline.h>
#include <gsl/gsl_linalg.h>

#define _DEBUG

namespace ub = boost::numeric::ublas;
using namespace std;

class CGForceMatching
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom)
    {
     n = SplineBond.GenerateGrid ( 0.08, 0.25, 0.1 ); // n - number of points
     n -= 1; // n - number of splines
  
     N = top->BeadCount();
     L = 0;        
     
        cout << "hey, someone wants to coarse grain\n";
        
        _A.resize(n-1, 2*(n+1), false); // to do: does it initialize _A with zeros?
        _b.resize(n-1, false);
        _A.clear();
        _b.clear();
        
   
        for (int i = 0; i < n-1; ++i) {
            _A(i,i) = SplineBond.A_prime(i, 0);
            _A(i, i+1) = SplineBond.B_prime(i,0) - SplineBond.A_prime(i,1);
            _A(i, i+2) = -SplineBond.B_prime(i,1);

            _A(i, n+1 + i) = SplineBond.C_prime(i,0);
            _A(i, n+1 + i+1) = SplineBond.D_prime(i,0) - SplineBond.C_prime(i,1);
            _A(i, n+1 + i+2) = -SplineBond.D_prime(i,1);
            
        }
    }
    
    void EndCG() {
        // Solving linear equations system

    //    double *pointer_m = new double [_b.size() * 2*(n+1)];
    //    double *pointer_b = new double [_b.size()];
        
        double* pointer_m = & _A(0,0);
        double* pointer_b = & _b(0);        
        
#ifdef _DEBUG        
        


        cout << _b.size() << " " << 2*(n+1) << endl;
        for(int i_line = 0; i_line < _b.size(); i_line++ ) {
            for(int j_col = 0; j_col < 2*(n+1); j_col++) {
//                cout << _A(i_line, j_col) << " ";
//                pointer_m[ i_line * 2*(n+1) + j_col ] = _A(i_line, j_col);
                cout << pointer_m[ i_line * 2*(n+1) + j_col ] << "\n";
            }
//            cout << "\n";
        }
        
        cout << _b.size() << endl;
        for(int i_line = 0; i_line < _b.size(); i_line++ ) {
//            cout << _b(i_line) << "\n";
//            pointer_b[ i_line ] = _b(i_line);
            cout << pointer_b[ i_line ] << endl;
        }        
        
#endif         


        
        gsl_matrix_view m
            = gsl_matrix_view_array (pointer_m, _b.size(), 2*(n+1));
    
        gsl_vector_view b
            = gsl_vector_view_array (pointer_b, _b.size());
    
        gsl_vector *x = gsl_vector_alloc (2*(n+1));

        gsl_vector *tau = gsl_vector_alloc (2*(n+1));

        gsl_vector *residual = gsl_vector_alloc (_b.size());
    
        gsl_linalg_QR_decomp (&m.matrix, tau);

        gsl_linalg_QR_lssolve (&m.matrix, tau, &b.vector, x, residual);

        /*
             
    gsl_matrix_view m 
      = gsl_matrix_view_array (_A, lin_dim, col_dim);

    gsl_vector_view b
      = gsl_vector_view_array (_b, lin_dim);
    
    gsl_vector *x = gsl_vector_alloc (col_dim);
    
    gsl_vector *tau = gsl_vector_alloc (col_dim);
        
    gsl_vector *residual = gsl_vector_alloc (lin_dim);
    
    gsl_linalg_QR_decomp (&m.matrix, tau);
    
    gsl_linalg_QR_lssolve (&m.matrix, tau, &b.vector, x, residual);
         
         */
        
        
        //output of the results 
        _x.resize(2*(n+1), false);
        _x.clear();
         
        cout << "write out results\n";
        for (int i =0 ; i < 2*(n+1); i++) {
            cout << "x[" << i <<"] = " << gsl_vector_get(x, i) << "\n";
            _x(i) = gsl_vector_get(x, i);
        }
         
        SplineBond.GetResult(& _x);
        SplineBond.PrintOutResult();
        
        gsl_vector_free (x);
        gsl_vector_free (tau);
        gsl_vector_free (residual);

    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
 //       cout << "yea, new configuration!\n";
        _A.resize(n-1 + 3*N*(L+1), 2*(n+1), true); // resize matrix
        _b.resize(n-1 + 3*N*(L+1), true);

        // here we have to set zero values to the new matrix elements - stupid thing
        for(int i_line = n-1 + 3*N*L; i_line < n-1 + 3*N*(L+1); i_line++ ) {
            for(int j_col = 0; j_col < 2*(n+1); j_col++) {
                _A(i_line, j_col) = 0;
            }
        }
        
#ifdef _DEBUG
      //  cout << "random matrix value: " << _A(n-1 + 3*N*L + 2*N, 2) << endl;
#endif             
        
        InteractionContainer &ic = conf->getTopology()->getBondedInteractions();
        InteractionContainer::iterator ia;

        // loop for the matrix:
        for(ia=ic.begin(); ia != ic.end(); ++ia) {
            if ((*ia)->BeadCount() == 2 ) {

               int i1 = (*ia)->getBeadId(0);
               int i2 = (*ia)->getBeadId(1);
               double r = (*ia)->EvaluateVar(*conf);
               int i = SplineBond.getInterval(r); 

               vec  n_ij = (*ia)->Grad(*conf, 0);

#ifdef _DEBUG
      //         cout << r << " " << -(conf->getF(i1) * n_ij) << endl;
#endif
               
                  _A(n-1 + 3*N*L+i1,i) = SplineBond.A(r)*n_ij.x(); 
                  _A(n-1 + 3*N*L+N+i1, i) = SplineBond.A(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i1, i) = SplineBond.A(r)*n_ij.z();  

                  _A(n-1 + 3*N*L+i1, i+1) = SplineBond.B(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i1, i+1) = SplineBond.B(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i1, i+1) = SplineBond.B(r)*n_ij.z();

                  _A(n-1 + 3*N*L+i1, n+i+1) = SplineBond.C(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i1, n+i+1) = SplineBond.C(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i1, n+i+1) = SplineBond.C(r)*n_ij.z();

                  _A(n-1 + 3*N*L+i1, n+1+i+1) = SplineBond.D(r)*n_ij.x(); 
                  _A(n-1 + 3*N*L+N+i1, n+1+i+1) = SplineBond.D(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i1, n+1+i+1) = SplineBond.D(r)*n_ij.z();


              n_ij = (*ia)->Grad(*conf, 1);

                  _A(n-1 + 3*N*L+i2,i) = SplineBond.A(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i2, i) = SplineBond.A(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i2, i) = SplineBond.A(r)*n_ij.z();

                  _A(n-1 + 3*N*L+i2, i+1) = SplineBond.B(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i2, i+1) = SplineBond.B(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i2, i+1) = SplineBond.B(r)*n_ij.z();

                  _A(n-1 + 3*N*L+i2, n+i+1) = SplineBond.C(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i2, n+i+1) = SplineBond.C(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i2, n+i+1) = SplineBond.C(r)*n_ij.z();

                  _A(n-1 + 3*N*L+i2, n+1+i+1) = SplineBond.D(r)*n_ij.x();
                  _A(n-1 + 3*N*L+N+i2, n+1+i+1) = SplineBond.D(r)*n_ij.y();
                  _A(n-1 + 3*N*L+2*N+i2, n+1+i+1) = SplineBond.D(r)*n_ij.z();


             }
//            cout << name << ": " << ((*ia)->EvaluateVar(*conf)) << endl;

        }

        // loop for the forces vector: 
        if ( conf->HasF() ) {
            vec Force(0., 0., 0.);
            for (int iatom = 0; iatom < N; ++iatom) {
                     Force = conf->getF(iatom);
                    _b(n-1 + 3*N*L + iatom) = Force.x();
                    _b(n-1 + 3*N*L + N+iatom) = Force.y();
                    _b(n-1 + 3*N*L + 2*N+iatom) = Force.z();
            }
        }
        else {
            cout << "ERROR: No forces in configuration!\n" ;   
        }
        L+=1; // update the frame counter
        
    }
    
protected:
    // _A*_x = _b
    
  ub::matrix<double> _A;
  ub::vector<double> _b; // F_ref
  ub::vector<double> _x; // 
    // _A(i, j) = 10;
    // _A.resize(n, m, true);
  int L; // counter for frames
  int N; //number of cg_beads
  int n; //number of splines
  CubicSpline SplineBond;
};


int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGForceMatching fmatch;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&fmatch);


    // initialize the readers/writers,
    // this will be combined in an initialize function later
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    
    // lets read in some program options
    namespace po = boost::program_options;
        
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    
    
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }
    // did an error occour?
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    return 0;
}

