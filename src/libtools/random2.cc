/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/tools/random2.h>
#include <stdexcept>
#include <iostream>
#include <votca/tools/constants.h>

namespace votca { namespace tools {

using namespace std;


void Random2::init( int nA1, int nA2, int nA3, int nB1 ) {

    nA1 = nA1 % 178 + 1;
    nA2 = nA2 % 178 + 1;
    nA3 = nA3 % 178 + 1;
    nB1 = nB1 % 169;

    if ((nA1 == 1) && (nA2 == 1) && (nA3 == 1)) {
        // Should not all be unity
        cout << flush
             << "WARNING: MARSAGLIA RNG INITIALISED INCORRECTLY. "
             << "ADAPTING SEEDS APPROPRIATELY."
             << endl;
        nA1 += nB1;
    }

    cout << flush
         << "INITIALIZED MARSFIELD WITH "
         << nA1 << " " << nA2 << " " << nA3 << " " << nB1
         << endl;

    /*
     initializes the global data of the
     MARSAGLIA pseudo random number generator
    */

    int mA1, mA2, mA3, mANEW, mB1, mHELP ;
    int i1, i2 ;
    double varS, varT ;
    MARSarray=std::vector<double>(MARS_FIELD_SIZE);
    mA1 = nA1 ; mA2 = nA2 ; mA3 = nA3 ;
    mB1 = nB1 ;
    MARSi  = 97 ; MARSj  = 33 ;

    for( i1 = 1 ; i1 < MARS_FIELD_SIZE ; i1++ )
    {
    varS = 0.0 ;
    varT = 0.5 ;
    for( i2 = 1 ; i2 < 25 ; i2++ )
    {
      mANEW = ((( mA1 * mA2 ) % 179 )*mA3 ) % 179 ;
      mA1 = mA2 ;
      mA2 = mA3 ;
      mA3 = mANEW ;
      mB1 = ( 53*mB1 + 1 ) % 169 ;
      mHELP = ( mB1 * mANEW ) % 64 ;
      if( mHELP > 31 )
        varS += varT ;
      varT *= 0.5 ;
    }

    MARSarray[ i1 ] = varS ;
    }

    MARSc  =   362436.0 / 16777216.0 ;
    MARScd =  7654321.0 / 16777216.0 ;
    MARScm = 16777213.0 / 16777216.0 ;

    return;
}

double Random2::rand_uniform( void ) {

    /*
     generates a pseudo random number 0 .. +1
     following the proposal of MARSAGLIA
    */

    double ranMARS ;

    ranMARS = MARSarray[ MARSi ] - MARSarray[ MARSj ] ;
    if( ranMARS < 0.0 )
    ranMARS += 1.0 ;

    MARSarray[ MARSi ] = ranMARS ;

    MARSi-- ;
    if( MARSi < 1 )
    MARSi = 97 ;

    MARSj-- ;
    if( MARSj < 1 )
    MARSj = 97 ;

    MARSc -= MARScd ;
    if( MARSc < 0.0 )
    MARSc += MARScm ;

    ranMARS -= MARSc ;
    if( ranMARS < 0.0 )
    ranMARS += 1.0 ;

    return ranMARS ;
}


/** generates a random integer number in the interval [0,max_int-1]
 */
int Random2::rand_uniform_int( int max_int ) {
    return floor (max_int * rand_uniform() );
}


/** generates a  gaussian distributed value 
 */
double Random2::rand_gaussian( double sigma ) {
  
    double r = sigma * sqrt( -2.0*log ( 1 - Random2::rand_uniform()) );
    double theta = 2.0 * conv::Pi * Random2::rand_uniform() ;
    return r * cos(theta); // second independent number is r*sin(theta)
}

}}
