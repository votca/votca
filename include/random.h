/*************************************************
     MARSAGLIA pseudo random number generator
     This function returns a double precision floating point number 
     uniformly distributed in the range [0,1)
*************************************************/

#ifndef _RANMARS_H_
#define _RANMARS_H_

#include <stdio.h>
#include <stdlib.h>

#define MARS_FIELD_SIZE 98

using namespace std;

class Random
{
public:
    static void init( int nA1, int nA2, int nA3, int nB1 );    
    static double rand( void );
    static void save( char *fileName );
    static void restore( char *fileName );

private:
    static double  *MARSarray, MARSc, MARScd, MARScm ;
    static int     MARSi, MARSj ;
 
};


#endif	/* _RANMARS_H_ */
