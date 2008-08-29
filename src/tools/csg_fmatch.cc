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
#include <stdio.h>

//#define _DEBUG

namespace ub = boost::numeric::ublas;
using namespace std;

class CGForceMatching
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom)
    {
        
        int lines_init = 0, colms_init = 0;  // initial size of _A 
        int sfnum; // number of spline functions for a cubicspline.        
        
        // set counters to zero value:
        line_cntr = col_cntr = 0;
        
        
        // SplineInfo for the first type of bond:
//        SplineInfo Bond1;
        Bond1.n = Bond1.Spline.GenerateGrid( 0.08, 0.37, 0.02) - 1;
        Bond1.bonded = true;
        Bond1.splineIndex = 0;
        Bond1.matr_pos = colms_init;
        
        Bond1.result.resize( 2*(Bond1.n + 1), false);
        Bond1.result.clear();
        
        //adjust initial matrix dimensions:
        lines_init += Bond1.n + 1;
        colms_init += 2 * (Bond1.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond1 );
        
        
        // SplineInfo for the second type of bond:
//        SplineInfo Bond2;
        Bond2.n = Bond2.Spline.GenerateGrid( 0.06, 0.42, 0.02) - 1;
        Bond2.bonded = true;
        Bond2.splineIndex = 1;
        Bond2.matr_pos = colms_init;
        
        Bond2.result.resize( 2*(Bond2.n + 1), false);
        Bond2.result.clear();
        
        //adjust initial matrix dimensions:
        lines_init += Bond2.n + 1;
        colms_init += 2 * (Bond2.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond2 );        
        
        // SplineInfo for the angle:
//        SplineInfo Angle;
        Angle.n = Angle.Spline.GenerateGrid( 113 * 0.0175, 126 * 0.0175, 2 * 0.0175 ) - 1;
        Angle.bonded = true;
        Angle.splineIndex = 2;
        Angle.matr_pos = colms_init;
        
        Angle.result.resize( 2*(Angle.n + 1), false);
        Angle.result.clear();
        
        //adjust initial matrix dimensions:
        lines_init += Angle.n + 1;
        colms_init += 2 * (Angle.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Angle );           
        
//        n = SplineBond.GenerateGrid ( 0.13, 0.17, 0.01 ); // n - number of points
//        n -= 1; // n - number of splines
  
        N = top->BeadCount();
        L = 0;        
     
        cout << "hey, someone wants to coarse grain\n";
        
        _A.resize(lines_init, colms_init, false); // to do: does it initialize _A with zeros?
        _b.resize(lines_init, false);
        _A.clear();
        _b.clear();        
        
        
//        _A.resize(n+1, 2*(n+1), false); // to do: does it initialize _A with zeros?
//        _b.resize(n+1, false);
//        _A.clear();
//        _b.clear();
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
        
            sfnum = (*is)->n;
            
            // smoothing conditions for first derivatives:
            // internal points
            
            for (int i =0; i < sfnum - 1; ++i) {
                _A(line_cntr + i, col_cntr + i) = (*is)->Spline.A_prime(i, 0);
                _A(line_cntr + i, col_cntr + i+1 ) = (*is)->Spline.B_prime(i,0) - (*is)->Spline.A_prime(i,1);
                _A(line_cntr + i, col_cntr + i+2 ) = -(*is)->Spline.B_prime(i,1);
                
                _A(line_cntr + i, col_cntr + sfnum+1 + i) = (*is)->Spline.C_prime(i,0);
                _A(line_cntr + i, col_cntr + sfnum+1 + i+1) = (*is)->Spline.D_prime(i,0) - (*is)->Spline.C_prime(i,1);
                _A(line_cntr + i, col_cntr + sfnum+1 + i+2) = -(*is)->Spline.D_prime(i,1);
            }
            
            // smoothing conditions for first derivatives:
            // external points 
            
            _A(line_cntr + sfnum-1, col_cntr + sfnum+1) = 1.0;
            _A(line_cntr + sfnum, col_cntr + 2*(sfnum+1) - 1) = 1.0;
            
            // update counters
            line_cntr += sfnum + 1;
            col_cntr += 2 * (sfnum + 1);
                
        }
        
/*        for (int i = 0; i < n-1; ++i) {
            _A(i,i) = SplineBond.A_prime(i, 0);
            _A(i, i+1) = SplineBond.B_prime(i,0) - SplineBond.A_prime(i,1);
            _A(i, i+2) = -SplineBond.B_prime(i,1);

            _A(i, n+1 + i) = SplineBond.C_prime(i,0);
            _A(i, n+1 + i+1) = SplineBond.D_prime(i,0) - SplineBond.C_prime(i,1);
            _A(i, n+1 + i+2) = -SplineBond.D_prime(i,1);
            
        }
        _A(n-1, n+1 ) = 1.0;
        _A(n, 2*(n+1) - 1) = 1.0;
*/        
    }
    
    void EndCG() {
        // Solving linear equations system
        
        // libgsl has problems with solving overdetermined systems by means of 
        // QR decomposition if the system doesn't have the full rank.
        // MATLAB can handle this problem easily.
        
        double* pointer_m = & _A(0,0);
        double* pointer_b = & _b(0);        
        
#ifdef _DEBUG        
        // this ifdef writes _A and _b with their dimensions to myfile.bin in binary format
        // this file can be used later to check the solution, e.g. with MATLAB
        // it also prints _A and _b to stdout.
        
        int lin_dim = _b.size();
        int col_dim = 2*(n+1);
        
        FILE * pFile;
        pFile = fopen ( "myfile.bin" , "wb" );
        fwrite(&lin_dim, sizeof(int), 1, pFile);
        fwrite(&col_dim, sizeof(int), 1, pFile);
        fwrite (pointer_m , sizeof(double) , 2*(n+1)*_b.size() , pFile );
        fwrite (pointer_b, sizeof(double), _b.size(), pFile);
        fclose (pFile);
        
        cout << _b.size() << " " << 2*(n+1) << endl;
        for(int i_line = 0; i_line < _b.size(); i_line++ ) {
            for(int j_col = 0; j_col < 2*(n+1); j_col++) {
                cout << pointer_m[ i_line * 2*(n+1) + j_col ] << "\n";
            }
            
        }
        
        cout << _b.size() << endl;
        for(int i_line = 0; i_line < _b.size(); i_line++ ) {
            cout << pointer_b[ i_line ] << endl;
        }        
        
#endif         

        gsl_matrix_view m
            = gsl_matrix_view_array (pointer_m, _b.size(), col_cntr);
    
        gsl_vector_view b
            = gsl_vector_view_array (pointer_b, _b.size());
    
        gsl_vector *x = gsl_vector_alloc (col_cntr);
        gsl_vector *tau = gsl_vector_alloc (col_cntr);       
        gsl_vector *residual = gsl_vector_alloc (_b.size());
    
        gsl_linalg_QR_decomp (&m.matrix, tau);

#ifdef _DEBUG 
        // prints TAU which is used by libgsl to stdout
        
        cout << "TAU ";
        for (int loop = 0; loop < 2*(n+1); loop++) {
             cout << gsl_vector_get(tau, loop) << " ";
        }
        cout << endl;
#endif
        
        
        gsl_linalg_QR_lssolve (&m.matrix, tau, &b.vector, x, residual);

        
#ifdef _DEBUG
        // prints RESIDUAL which is used by libgsl to stdout
        cout << "RESIDUAL \n";
        for (int loop = 0; loop < _b.size(); loop++) {
             cout << gsl_vector_get(residual, loop) << endl;
        } 
#endif
        
        //output of the results 
        _x.resize(col_cntr, false);
        _x.clear();
         
        cout << "write out results\n";
        // do we really need this _x vector?
        // we need many of them for each spline instead.
        // they are implemented in SplineInfo
/*        for (int i =0 ; i < col_cntr; i++) {
            cout << "x[" << i <<"] = " << gsl_vector_get(x, i) << "\n";
            _x(i) = gsl_vector_get(x, i);
        }
*/        
        cout << "r  " << "F(r)  " << endl;
        
        // have to change it - it's bullshit
//        Bond1.Spline.GetResult(& _x);
//        Bond1.Spline.PrintOutResult();
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
            int &mp = (*is)->matr_pos;
            int &nsf = (*is)->n;
            cout << "interaction No. " << (*is)->splineIndex << endl;
            
            for (int i = 0; i < 2*(nsf + 1); i++ ) {
               (*is)->result[i] = gsl_vector_get(x, i + mp);
            }
            (*is)->Spline.GetResult( & (*is)->result );
            (*is)->Spline.PrintOutResult();
        }        
        
        gsl_vector_free (x);
        gsl_vector_free (tau);
        gsl_vector_free (residual);
 
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
 //       cout << "yea, new configuration!\n";
        _A.resize(line_cntr + 3*N*(L+1), col_cntr, true); // resize matrix _A
        _b.resize(line_cntr + 3*N*(L+1), true);          // resize vector _b

//        _A.resize(n+1 + 3*N*(L+1), 2*(n+1), true); // resize matrix _A
//        _b.resize(n+1 + 3*N*(L+1), true);          // resize vector _b

        // here we have to set zero values to the new matrix elements - stupid thing
        for(int i_line = line_cntr + 3*N*L; i_line < line_cntr + 3*N*(L+1); i_line++ ) {
            for(int j_col = 0; j_col < col_cntr; j_col++) {
                _A(i_line, j_col) = 0;
            }
        }
                  
        InteractionContainer &ic = conf->getTopology()->getBondedInteractions();
        InteractionContainer::iterator ia;

        // loop for the matrix:
        for(ia=ic.begin(); ia != ic.end(); ++ia) {
            
               int beads_in_int = (*ia)->BeadCount(); // 2 for bonds, 3 for angles, 4 for dihedrals
                
               int index = (*ia)->getGroupId(); // unique for every interaction type

               CubicSpline &SP = Splines[ index ]->Spline;
               int  &mpos = Splines[ index ]->matr_pos;
               int  &nsp = Splines[ index ]->n;
               
               double var = (*ia)->EvaluateVar(*conf); // value of bond, angle, or dihedral
               int i = SP.getInterval(var);   // corresponding spline interval

               for (int loop = 0; loop < beads_in_int; loop ++) {
                   int ii = (*ia)->getBeadId(loop);
                   vec gradient = (*ia)->Grad(*conf, loop);
                   
                  _A(line_cntr + 3*N*L+ii, mpos + i) += SP.A(var)*gradient.x(); 
                  _A(line_cntr + 3*N*L+N+ii, mpos + i) += SP.A(var)*gradient.y();
                  _A(line_cntr + 3*N*L+2*N+ii, mpos + i) += SP.A(var)*gradient.z();  

                  _A(line_cntr + 3*N*L+ii, mpos + i+1) += SP.B(var)*gradient.x();
                  _A(line_cntr + 3*N*L+N+ii, mpos + i+1) += SP.B(var)*gradient.y();
                  _A(line_cntr + 3*N*L+2*N+ii, mpos + i+1) += SP.B(var)*gradient.z();

                  _A(line_cntr + 3*N*L+ii, mpos + nsp+i+1) += SP.C(var)*gradient.x();
                  _A(line_cntr + 3*N*L+N+ii, mpos + nsp+i+1) += SP.C(var)*gradient.y();
                  _A(line_cntr + 3*N*L+2*N+ii, mpos + nsp+i+1) += SP.C(var)*gradient.z();

                  _A(line_cntr + 3*N*L+ii, mpos + nsp+1+i+1) += SP.D(var)*gradient.x(); 
                  _A(line_cntr + 3*N*L+N+ii, mpos + nsp+1+i+1) += SP.D(var)*gradient.y();
                  _A(line_cntr + 3*N*L+2*N+ii, mpos + nsp+1+i+1) += SP.D(var)*gradient.z();                   
                   
               }
        }

        // loop for the forces vector: 
        if ( conf->HasF() ) {
            vec Force(0., 0., 0.);
            for (int iatom = 0; iatom < N; ++iatom) {
                     Force = conf->getF(iatom);
                    _b(line_cntr + 3*N*L + iatom) = Force.x();
                    _b(line_cntr + 3*N*L + N+iatom) = Force.y();
                    _b(line_cntr + 3*N*L + 2*N+iatom) = Force.z();
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
  
  int line_cntr, col_cntr; // counters for lines and coloumns in _A 
  
//  int n; //number of splines
//  CubicSpline SplineBond;
  
  struct SplineInfo {
        int n; //number of splines
        int splineIndex; // interaction index for bonded interactions
        bool bonded;     // true for bonded interactions, false for non-bonded
        CubicSpline Spline;
        int matr_pos;    // position in the _A matrix (first coloumn which is occupied with
                         // this particular spline
        ub::vector<double> result;
  };
  SplineInfo Bond1;
  SplineInfo Bond2;
  SplineInfo Angle;
  
  typedef vector<SplineInfo *> SplineContainer;
  SplineContainer Splines;
  
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

