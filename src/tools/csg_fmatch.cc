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
#include <sstream>

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
        const int N_frames = 50001; // Number of frames in the trajectory
        
        // set counters to zero value:
        line_cntr = col_cntr = 0;
        
        
        // SplineInfo for the first type of bond:
        Bond1.n = Bond1.Spline.GenerateGrid(  0.256, 0.337, 0.005 ) - 1;
        Bond1.bonded = true;
        Bond1.splineIndex = 0;
        Bond1.splineName = "bond1";
        Bond1.matr_pos = colms_init;
        
        Bond1.result.resize( 2*(Bond1.n + 1), false);
        Bond1.result.clear();
        
        //adjust initial matrix dimensions:
        lines_init += Bond1.n + 1;
        colms_init += 2 * (Bond1.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond1 );
        
        
        // SplineInfo for the second type of bond:
/*        Bond2.n = Bond2.Spline.GenerateGrid( 0.345, 0.395, 0.01) - 1;
        Bond2.bonded = true;
        Bond2.splineIndex = 1;
        Bond2.splineName = "bond2";
        Bond2.matr_pos = colms_init;
        
        Bond2.result.resize( 2*(Bond2.n + 1), false);
        Bond2.result.clear();
        
        //adjust initial matrix dimensions:
        lines_init += Bond2.n + 1;
        colms_init += 2 * (Bond2.n + 1);
                
        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Bond2 );        
*/        
        // SplineInfo for the angle:
        Angle1.n = Angle1.Spline.GenerateGrid( 69.8 * 0.0175, 170 * 0.0175, 5 * 0.0175 ) - 1;
        Angle1.bonded = true;
        Angle1.splineIndex = 1;
        Angle1.splineName = "angle1";
        Angle1.matr_pos = colms_init;

        Angle1.result.resize( 2*(Angle1.n + 1), false);
        Angle1.result.clear();

        //adjust initial matrix dimensions:
        lines_init += Angle1.n + 1;
        colms_init += 2 * (Angle1.n + 1);

        //Add SplineInfo to SplineContainer:
        Splines.push_back( &Angle1 );
        
        // angle 2
        
            
//        n = SplineBond.GenerateGrid ( 0.13, 0.17, 0.01 ); // n - number of points
//        n -= 1; // n - number of splines
  
        N = top->BeadCount(); // Number of beads in topology
        L = 0;                // Initial frame in trajectory  
        cout << "hey, someone wants to coarse grain\n";
        
        
        // B_constr matrix contains continuity conditions for the spline first
        // derivatives.
        B_constr.resize(lines_init, colms_init, false);
        B_constr.clear();
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
        
            sfnum = (*is)->n;            
            (*is)->Spline.AddBCToFitMatrix(B_constr, line_cntr, col_cntr);
            
            // update counters
            line_cntr += sfnum + 1;
            col_cntr += 2 * (sfnum + 1);
                
        }
        
        _A.resize( 3*N*N_frames, col_cntr, false); // resize matrix _A
        _b.resize( 3*N*N_frames, false);          // resize vector _b   
        _A.clear();
        _b.clear();        
               
    }
    
    void EndCG() {
        string force_raw = "_force_raw.dat";
        char file_name[20];
        double accuracy; // accuracy for output. Should be different for bonds and angles.
        
        ofstream out_file;
         
        // Solving linear equations system
        
        _x.resize(col_cntr);
        _x.clear();

        ub::matrix<double> Q; 
        Q.resize(col_cntr, col_cntr, false );
        Q.clear();
        
        ub::matrix<double> A2;  
        A2.resize(_A.size1(), col_cntr / 2, false );
        A2.clear();
                
        ub::matrix<double> Q_k;
        Q_k.resize(col_cntr, col_cntr, false);
        Q_k.clear();
       
        ub::identity_matrix<double> I (col_cntr); 

        ub::vector<double> v; 
        v.resize(col_cntr, false);
        v.clear();        
        
        // To proceed we need to factorize B^T = Q*R. We need matrix Q for further
        // calculations
        B_constr = trans(B_constr);
        
        double* pointer_Bcnstr = & B_constr(0,0); 
        
        gsl_matrix_view B_t 
          = gsl_matrix_view_array (pointer_Bcnstr, col_cntr, line_cntr);
     
        gsl_vector *tau = gsl_vector_alloc (line_cntr); 
    
        gsl_linalg_QR_decomp (&B_t.matrix, tau);   

        // Extraction of Q matrix from tau and B_t, where it is stored in a tricky way.
        Q = I;
              
        for (int k = line_cntr; k > 0 ; k--) {
           
            for (int icout = 0; icout < k - 1; icout++) {
                 v(icout) = 0;
            }
            v(k - 1) = 1.0;

            for (int icout = k; icout < col_cntr; icout++) {
                 v(icout) = gsl_matrix_get(&B_t.matrix, icout, k - 1 );
            }
           
            Q_k = I - gsl_vector_get(tau, k - 1 ) * outer_prod ( v, v );
            Q = prec_prod(Q, Q_k);
           
        }
    
            
        Q = trans(Q);     
    
        // Calculate _A * Q and store the result in _A
        _A = prec_prod(_A, Q);
    
        // _A = [A1 A2], so A2 is just a block of _A
        for (int iraw = 0; iraw < _A.size1(); iraw++) {
             for (int icol = _A.size2() / 2; icol < _A.size2(); icol++) {
                A2(iraw, icol - _A.size2() / 2) = _A(iraw, icol);
             }
        }    
        
        
  
        double* pointer_m = & A2(0,0);
        double* pointer_b = & _b(0);        
        
        gsl_matrix_view m
            = gsl_matrix_view_array (pointer_m, A2.size1(), A2.size2() );
    
        gsl_vector_view b
            = gsl_vector_view_array (pointer_b, A2.size1());
    
        gsl_vector *x = gsl_vector_alloc ( A2.size2() );
        gsl_vector *tau2 = gsl_vector_alloc ( A2.size2() );       
        gsl_vector *residual = gsl_vector_alloc ( A2.size1() );
    
        gsl_linalg_QR_decomp (&m.matrix, tau2);
        
        gsl_linalg_QR_lssolve (&m.matrix, tau2, &b.vector, x, residual);
        
        for (int i = 0; i < col_cntr / 2; i++ ) {
               _x[i] = 0.0;
        }  
    
        for (int i = col_cntr / 2; i < col_cntr; i++ ) {
               _x[i] = gsl_vector_get(x, i - col_cntr / 2 );
        }      
    
        // To get the final answer this vector should be multiplied by matrix Q
        _x = prec_prod( Q, _x );        
        
        
        SplineContainer::iterator is;
                
        for(is=Splines.begin(); is != Splines.end(); ++is) {
            int &mp = (*is)->matr_pos;
            int &nsf = (*is)->n;
            
            file_name[0] = '\0';
            strcpy(file_name, ((*is)->splineName).c_str() );
            strcat(file_name, force_raw.c_str());
            out_file.open(file_name);
            
            out_file << "# interaction No. " << (*is)->splineIndex << endl;
            
            for (int i = 0; i < 2*(nsf + 1); i++ ) {
//               (*is)->result[i] = gsl_vector_get(x, i + mp);
                (*is)->result[i] = _x[ i + mp ];
            }
            //(*is)->Spline.GetResult( & (*is)->result );
            //(*is)->Spline.PrintOutResult();
            (*is)->Spline.setSplineData( (*is)->result );
            
            if ( ((*is)->splineName)[0] == 'a' ) accuracy = 0.05;
            else if ( ((*is)->splineName)[0] == 'b' ) accuracy = 0.001;
            
            (*is)->Spline.Print(out_file, accuracy);
            
            out_file.close();
        }        
        
        gsl_vector_free (x);
        gsl_vector_free (tau);
        gsl_vector_free (residual);
 
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
                  
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
                  
                   SP.AddToFitMatrix(_A, var, 
                           3*N*L + ii, mpos, gradient.x());
                   SP.AddToFitMatrix(_A, var, 
                           3*N*L + N + ii, mpos, gradient.y());
                   SP.AddToFitMatrix(_A, var,
                           3*N*L + 2*N + ii, mpos, gradient.z());                
               }
        }

        // loop for the forces vector: 
        if ( conf->HasF() ) {
            vec Force(0., 0., 0.);
            for (int iatom = 0; iatom < N; ++iatom) {
                     Force = conf->getF(iatom);
                    _b( 3*N*L + iatom) = Force.x();
                    _b( 3*N*L + N+iatom) = Force.y();
                    _b( 3*N*L + 2*N+iatom) = Force.z();
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
  ub::matrix<double> B_constr;
    // _A(i, j) = 10;
    // _A.resize(n, m, true);
  int L; // counter for frames
  int N; //number of cg_beads
  
  int line_cntr, col_cntr; // counters for lines and coloumns in B_constr 
  
  struct SplineInfo {
        int n; //number of splines
        int splineIndex; // interaction index for bonded interactions
        bool bonded;     // true for bonded interactions, false for non-bonded
        CubicSpline Spline;
        int matr_pos;    // position in the _A matrix (first coloumn which is occupied with
                         // this particular spline
        ub::vector<double> result;
        string splineName;
  };
  SplineInfo Bond1;
//  SplineInfo Bond2;
  SplineInfo Angle1;
//  SplineInfo Angle2;
//  SplineInfo Angle3;
  
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

