/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#include "random.h"

double  *Random::MARSarray, Random::MARSc, Random::MARScd, Random::MARScm ;
int     Random::MARSi, Random::MARSj ;

void Random::init( int nA1, int nA2, int nA3, int nB1 )
{
  /*
     initializes the global data of the
     MARSAGLIA pseudo random number generator
  */

  int mA1, mA2, mA3, mANEW, mB1, mHELP ;
  int i1, i2 ;
  double varS, varT ;

  MARSarray = (double*)malloc(MARS_FIELD_SIZE*sizeof(double)); 

  mA1 = nA1 ; mA2 = nA2 ; mA3 = nA3 ;
  mB1 = nB1 ;
  MARSi  = 97 ; MARSj  = 33 ;

  for( i1 = 1 ; i1 < 98 ; i1++ )
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

void Random::save( char *fileName )
{
  	FILE *ranFP;
  	int c[2];
  	double w[3];

  	ranFP = fopen(fileName, "wb");
  	if (ranFP==NULL)
  	{
		printf("Dateifehler in 'restore_RAN'\n");
		exit(-1);
	}
  
	  c[0] = MARSi; c[1] = MARSj;
  	fwrite(c, sizeof(int), 2, ranFP);
  	w[0] = MARSc; w[1] = MARScd; w[2] = MARScm;
  	fwrite(w, sizeof(double), 3, ranFP);
  	fwrite(MARSarray, sizeof(double), MARS_FIELD_SIZE, ranFP);
  	fclose(ranFP);
}

void Random::restore( char *fileName )
{
	FILE *ranFP;
  	double w[3];
  	int c[2];

  	ranFP = fopen(fileName, "rb");
  	if (ranFP==NULL)
  	{
		printf("Dateifehler in 'restore_RAN'\n");
		exit(-1);
	}
  
  	  fread(c, sizeof(int), 2, ranFP);
    	MARSi = c[0]; MARSj = c[1];
    	fread(w, sizeof(double), 3, ranFP);
    	MARSc = w[0]; MARScd = w[1]; MARScm = w[2];
    	fread(MARSarray, sizeof(double), MARS_FIELD_SIZE, ranFP);
    	fclose(ranFP);
}

double Random::rand_uniform( void )
{
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
int Random::rand_uniform_int( int max_int ) 
{
  return floor (max_int * rand_uniform() );
}

/** generates a  gaussian distributed value 
 */
double Random::rand_gaussian( double sigma )
{
    double r = sigma * sqrt( -2.0*log ( 1 - Random::rand_uniform()) );
    double theta = 2.0 * _pi * Random::rand_uniform() ;
    return r * cos(theta); // second independent number is r*sin(theta)
}
