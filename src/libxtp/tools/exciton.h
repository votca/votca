/* 
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_XTP_EXCITON_H
#define _VOTCA_XTP_EXCITON_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/qmtool.h>
#include <votca/ctp/segment.h>
#include <votca/tools/constants.h>

#include <votca/tools/linalg.h>
#include <votca/tools/constants.h>
namespace votca { namespace xtp {
    using namespace std;
    
class Exciton : public ctp::QMTool
{
public:

    Exciton() { };
   ~Exciton() { };

    string Identify() { return "exciton"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    
 

private:
    
    string      _orbfile; // file containining the MOs from qmpackage...
    string      _logfile; // file containining the Energies etc... from qmpackage...
    string      _xyzfile;

    string      _package;
    Property    _package_options;
    Property    _gwbse_options;
    Property    _summary;
    
    string      _output_file; // .orb file to parse to
    string      _xml_output;    // .xml output

    string      _reporting;    
    ctp::Logger      _log;
    
    bool _do_dft_input;
    bool _do_dft_run;
    bool _do_dft_parse;
    bool _do_gwbse;
    bool _do_optimize;
    
    int _opt_state;
    double _displacement;
    double _convergence;
    double _trust_radius;
    double _trust_radius_max;
    double _delta_energy_estimate;
    double _norm_delta_pos;
    string _spintype;
    string _forces; 
    string _opt_type;    
    
    string _guess_orbA;
    string _guess_orbB;
    bool _do_guess;
    
    int _natoms;
    int _iteration;
    ub::matrix<double> _force;
    ub::matrix<double> _force_old;
    ub::matrix<double> _xyz_shift;
    ub::matrix<double> _speed;
    ub::matrix<double> _current_xyz;
    ub::matrix<double> _old_xyz; 
    ub::matrix<double> _trial_xyz; 
    ub::matrix<double> _hessian;
    
    bool _step_accepted;
    bool _update_hessian;
    bool _restart_opt;
    
    void ExcitationEnergies( QMPackage* _qmpackage, vector <ctp::Segment* > _segments, Orbitals* _orbitals );
    void ReadXYZ( ctp::Segment* _segment, string filename);
    void Orbitals2Segment(ctp::Segment* _segment, Orbitals* _orbitals);
    void Coord2Segment(ctp::Segment* _segment );
    

    void PrepareGuess(         Orbitals *_orbitalsA, 
                               Orbitals *_orbitalsB, 
                               Orbitals *_orbitalsAB);

    void OrthonormalizeGuess ( ctp::Segment* _segment, Orbitals* _orbitals );
    
    
    void BFGSStep( int& _iteration, bool& _update_hessian,  ub::matrix<double>& _force, ub::matrix<double>& _force_old,  ub::matrix<double>& _current_xyz, ub::matrix<double>&  _old_xyz, ub::matrix<double>& _hessian ,ub::matrix<double>& _xyz_shift ,ub::matrix<double>& _trial_xyz  );
    void ReloadState();
    void NumForceForward(double energy, vector <ctp::Atom* > _atoms, ub::matrix<double>& _force, QMPackage* _qmpackage,vector <ctp::Segment* > _segments, Orbitals* _orbitals );
    void NumForceCentral(double energy, vector <ctp::Atom* > _atoms, ub::matrix<double>& _force, QMPackage* _qmpackage,vector <ctp::Segment* > _segments, Orbitals* _orbitals );
    
    void WriteIteration( FILE* out, int _iteration, ctp::Segment* _segment, ub::matrix<double>& _force  );
    
    string Convergence( bool _converged ) { 
        
        string _converged_string;
        if ( _converged )  _converged_string = " (converged)";
        if ( !_converged )  _converged_string = " (not converged)";
        
        
        return _converged_string;
    }


};



}}


#endif
