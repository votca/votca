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

double Random::rand( void )
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
