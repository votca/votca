/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#ifndef _CALC_INTEGRALS_DFT_H
#define	_CALC_INTEGRALS_DFT_H

#include <votca/ctp/parallelpaircalc.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/gaussian.h>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <sys/stat.h>

namespace votca { namespace ctp {
    
/**
* \brief Density-functional-based electronic coupling elements
*
* DFT-based electronic coupling elements for all conjugated
* segments from the neighbor list. Requires molecular orbitals in GAUSSIAN
* or TURBOMOLE format.
*
* Callname: idft
*/

class IDFT : public ParallelPairCalculator
{
public:

    IDFT() {};
   ~IDFT() {};
   
    // Topology is not needed; shall removed later with a different Application class
    void    Initialize(ctp::Topology *top, tools::Property *options );
    
    string  Identify() { return "IDFT"; }
    
    double getCouplingElement( int level1, int level2 );   
    void EvalPair(Topology *top, QMPair *qmpair, int slot);
    
/*  
    void    EvalPair(Topology *top, QMPair *pair, int slot);
    void    CleanUp();
*/

private:

    static const double _conv_Hrt_eV = 27.21138386;

    bool                _has_integrals;
    bool                _has_degeneracy;
    
    int                 _max_occupied_levels;
    int                 _max_unoccupied_levels;

    string              _package;
    Property            _package_options; 
    
    double              _energy_difference;    
    
    string              _orbitalsA_file;
    string              _orbitalsB_file;
    string              _orbitalsAB_file;
    
    string              _outParent;
        
    Orbitals            _orbitalsA;
    Orbitals            _orbitalsB;
    Orbitals            _orbitalsAB;

    ub::matrix<double>  _JAB;     
    
    void        SQRTOverlap(ub::symmetric_matrix<double> &S, ub::matrix<double> &Sm2);
    void        CalculateIntegrals();   
    void        ParseOptionsXML( tools::Property *opt);    
    
};

}}
#endif	/* _CALC_INTEGRALS_DFT_H */
