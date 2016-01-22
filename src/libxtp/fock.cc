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

#include <votca/xtp/fock.h>
#include <votca/tools/property.h>

namespace votca { namespace xtp {

//ZINDO/S parameters for first three rows
//Slater exponenets are assumed to be the same for s/p.  This is only the case up to row 3.  Including transition metals would require copious
//rewriting.
//jjk: Add Zindo_s parameters for Al (from Orca) - Values for H,C,N,O all agree.  Haven't yet checked for other atoms

/* data now read in from file
static double Beta_zindo_s[18] = {-12.000000, 0.000000, 0.000000, 0.000000, 0.000000, -17.000000, -26.000000, -34.000000, -44.000000, 0.000000, 0.000000, -6.000000, -11.3000000, 0.000000, 0.000000, -15.000000, -19.000000, 0.00000 } ; // JK : made these static so that only 1 copy exists for all instantiated class objects

static double Mu_zindo_s  [18]  ={ 1.200000, 0.000000, 0.000000, 0.000000, 0.000000, 1.625000, 1.950000, 2.275000, 2.600000, 0.000000, 0.000000, 1.103000, 1.370000, 0.000000, 0.000000, 1.816670, 2.033330, 0.000000 };   //Zindo_s parameters for the first three rows

static double Beta_zindo_1[18]  ={-11.800000,-0.000000,-9.000000,-13.000000,-17.000000,-21.000000,-25.000000,-31.000000,-42.000000,-0.000000,-7.720000,-9.450000,-11.300000,-13.000000,-15.100000,-16.000000,-18.000000,-0.000000  };

static double Mu_zindo_1  [18]  ={1.200000,1.700000,0.650000,0.975000,1.300000,1.625000,1.950000,2.275000,2.600000,2.925000,0.733330,1.103000,1.166670,1.383330,1.600000,1.816670,2.033330,2.250000  };  //Zindo_1  --''--

*/

 					          //!!N.B.!! :
					          //Slater exponenets are assumed to be the same for s/p.  
						  //This is only the case up to row 3.
						  //Including transition metals would require a major re-write.
const static int HASH_S = 100000; // resolution in hash table to STO overlap of AOs
static double _1S1S [HASH_S+1];
static double _2S2S [HASH_S+1];
static double _2S2Sig [HASH_S+1];
static double _2Sig2Sig [HASH_S+1];
static double _2Pi2Pi [HASH_S+1];
static double _1S2S_hc [HASH_S+1];
static double _1S2Sig_hc [HASH_S +1];


const static double PMAX=1000.0; // this defines the maximum p - refer to Mulliken paper for more detail
const static double PMIN=0.0;

bool fock::_init_HASH=false;

/*
void fock::set_zindo_1() {
	Beta=Beta_zindo_1;
	Mu=Mu_zindo_1;
    }
void fock::set_zindo_s() {
	Beta=Beta_zindo_s;
	Mu=Mu_zindo_s;
    }
         * */
void fock::SetParameters(){
    char * votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL)
        throw std::runtime_error("VOTCASHARE not set, cannot open INDO parameters.");
    string nameParameter = string(getenv("VOTCASHARE"))+string("/xml/xtp/INDOParameters.xml");
    Property Options;
    load_property_from_xml(Options,nameParameter);
    
    Property values = (Options.get("INDOparam"));
    
    string muS, betaS;
    muS    = values.get("mu").as<string>();
    betaS  = values.get("beta").as<string>();
    
    Tokenizer tokM(muS, "\n\t ");
    Tokenizer tokB(betaS, "\n\t ");
    vector<string> Mutok, Betatok;
    tokM.ToVector(Mutok);
    tokB.ToVector(Betatok);
    vector<string>::iterator it= Mutok.begin();
    for (; it != Mutok.end();++it)
    {
        Mu.push_back(boost::lexical_cast<double>(*it));
    }
     it= Betatok.begin();
    for (; it != Betatok.end();++it)
    {
        Beta.push_back(boost::lexical_cast<double>(*it));
    }
}

void fock::CheckParameters(const mol_and_orb & A ){
    int Nat = A.getN();
    for (int i=0;i<Nat;i++){
        unsigned int lbl = A.getlbl(i);
        if (lbl>=Mu.size() ){
            throw runtime_error("not enough specifications for mu ");
        }
        if (lbl>=Beta.size() ){
            throw runtime_error("not enough specifications for beta ");
        }
        if (Mu[lbl] ==0){
            throw runtime_error("atom with lbl "+boost::lexical_cast<string>(lbl)+ " has no mu param");
        }
        if (Beta[lbl] ==0){
            throw runtime_error("atom with lbl "+boost::lexical_cast<string>(lbl)+ " has no beta param");
        }
    }
    
}

inline double Abs (const double & a)
{
    if (a< 0.0) return -a;
    return a;
}
inline double A6(double *p_p, const double & EXP){
    	return EXP *(p_p[0] + 6.0*p_p[1] + 30.0*p_p[2] + 120.0 *p_p[3] + 360.0*p_p[4] + 720.0*p_p[5] + 720*p_p[6]);
}
	
inline double A5( double *p_p, const double & EXP){
        return EXP*(p_p[0] + 5.0*p_p[1] + 20.0 *p_p[2] + 60*p_p[3] + 120.0 * p_p[4] + 120.0 *p_p[5]);
}
inline double A4( double * p_p, const double & EXP)
{
	return EXP*(p_p[0] + 4.0*p_p[1] + 12.0*p_p[2] + 24.0*p_p[3] + 24.0*p_p[4] );
}
inline double A3(double * p_p, const double & EXP)
{
	return EXP*( p_p[0] + 3.0 * p_p[1] + p_p[2] *6.0 + p_p[3]*6.0);
}
inline double A2(double * p_p, const double & EXP)
{
	return EXP*( 2.0*p_p[2] + 2.0*p_p[1] + p_p[0]);
}

inline double A1(double * p_p, const double & EXP)
{
	return EXP*(p_p[0]+p_p[1]);
}
inline double A0(double * p_p, const double & EXP)
{
//	cout << EXP << " " << p_p[0]
	return EXP*p_p[0];
}

inline double B6(const double & cos_h, const double &sin_h, double *p_pt){
    return sin_h *(2.0*p_pt[0] + 60.0 *p_pt[2] + 720.0 * p_pt[4] +1440.0 * p_pt[6]) - cos_h * (12.0*p_pt[1] + 240.0 *p_pt[3] +1440.0 * p_pt[5]  );
}

inline double B5(const double & cos_h, const double &sin_h, double *p_pt ){
    return sin_h * ( 10.0*p_pt[1] + 120.0 *p_pt[3]+240.0*p_pt[5]  ) - cos_h * (2.0*p_pt[0] + 40*p_pt[2] + 240.0 * p_pt[4]);
}

inline double B4(const double & cos_h, const double & sin_h, double *p_pt)
{
	return sin_h * (2.0*p_pt[0] + 24.0*p_pt[2] + 48.0*p_pt[4] ) - cos_h * (8.0*p_pt[1] + 48.0*p_pt[3]) ;
}
inline double B3(const double & cos_h,const double & sin_h, double *p_pt)
{
	return sin_h * ( 6.0*p_pt[1] + 12.0 *p_pt[3] ) - cos_h * (2.0*p_pt[0] + 12.0 *p_pt[2]);
}
inline double B2( const double & cos_h,const double & sin_h, double *p_pt)
{
//	cout << sin_h*(2.0*p_pt[0]+4.0*p_pt[2]) << " " << -cos_h*4.0*p_pt[1] << endl;
	return sin_h*(2.0*p_pt[0]+4.0*p_pt[2])-cos_h*4.0*p_pt[1];
}
inline double B1(const double & cos_h, const double & sin_h, double *p_pt)
{
//	cout << sin_h*2.0 << " " << p
	return sin_h*2.0*p_pt[1] - 2.0*cos_h*p_pt[0];
}

inline double B0(const double &cos_h, const double & sin_h, double *p_pt){
    	return sin_h*2.0*p_pt[0];
}
inline double B0(const double & sin_h, double *p_pt)
{
//	cout << sin_h * 2.0 << " " << p_pt[0] << endl; //checked
	return sin_h*2.0*p_pt[0];
}
/////////////////////////////////////////////////////////////////////////
//

int fock::init_HASH(){
    double p;
    double t, pt;
    double p_int=(PMAX-PMIN)/HASH_S;
    int i=0;
    double p_p[5];
    double p_pt[5];
    double EXP, cos_h, sin_h;
    p=p_int/2.0;
    t=(-0.425/2.825);
    for ( i=0; i< HASH_S; i++)
    {
	pt =p*t;
	EXP = exp( -p);
	p_p[0]=1.0/p;
	p_pt[0]=1.0/pt;
	for(int j=1;j<5;j++) p_pt[j] = p_pt[j-1] /pt;
	for(int j=1;j<5;j++) p_p[j] = p_p[j-1] / p;
	cos_h = cosh(pt);
	sin_h = sinh(pt);
	_1S1S[i]     = p*p*p / 6.0 * (A2(p_p, EXP) *3.0 - A0(p_p, EXP));
	_2S2S[i]     = EXP * (1.0 + p + 4.0 *p*p / 9.0 + p*p*p / 9.0 + p*p*p*p/45.0);
	_2S2Sig[i]   = p*p*p*p*p / ( 60.0 * sqrt(3.0) ) * ( 5*A3(p_p, EXP) - A1( p_p, EXP) );
        _2Sig2Sig[i] = p*p*p*p*p / 120.0 * ( 5.0 * A4(p_p, EXP) - 18.0 * A2(p_p, EXP) + 5.0 * A0(p_p, EXP) );
	_2Pi2Pi[i]   =  p*p*p*p*p / 120.0 * ( 5.0 * A4(p_p, EXP) -  6.0 * A2(p_p, EXP) + A0(p_p, EXP) );
	_1S2S_hc[i]  =  p*p*p*p / ( 8.0 * sqrt(3.0 )) * pow(1-t, 2.5) * pow(1+t, 1.5) *
             ( A3(p_p,EXP) * B0(sin_h, p_pt)        - A2(p_p, EXP) * B1(cos_h, sin_h, p_pt) -
               A1(p_p,EXP) * B2(cos_h, sin_h, p_pt) + A0(p_p, EXP) * B3(cos_h, sin_h, p_pt) )  ;
	_1S2Sig_hc[i] = p*p*p*p / 8.0 * pow(1-t, 2.5) * pow(1+t, 1.5) *
          (-A3(p_p,EXP) * B1(cos_h, sin_h, p_pt) + A2(p_p, EXP) * B0(sin_h, p_pt)       +
           A1(p_p,EXP) * B3(cos_h, sin_h, p_pt) - A0(p_p, EXP) * B2(cos_h, sin_h, p_pt) ) ;
	p+=p_int;
    }
    _1S1S[HASH_S]=0.0;
    _2S2S[HASH_S]=0.0;
    _2S2Sig[HASH_S]=0.0;
    _2Sig2Sig[HASH_S]=0.0;
    _2Pi2Pi[HASH_S]=0.0;
    _1S2S_hc[HASH_S]=0.0;
    _1S2Sig_hc[HASH_S]=0.0;
    return 0;
}

void fock::calc_F_el_with_HASH( const double & p,  double & t,  const double & beta, const double & x, const double & y, const double & z, const int & c1, const int & c2, const int &lbl1, const int &lbl2) const  {

static int i;
static int P;
if ( p > PMAX) P=HASH_S;
else P=int ( ((p-PMIN)/(PMAX-PMIN))*double (HASH_S) );
//these are the sto overlap factors//
static double IsIs, IsIIsig,IsIIs, IsIIIsig, IsIIIs;
static double IIsIIs, IIsIIsig, IIsigIIs, IIsIIIs, IIsIIIsig;
static double IIsigIIsig,IIsigIIIsig, IIpiIIpi, IIpiIIIpi, IIsigIIIs, IIIsigIIs, IIIsIIsig;
static double IIIsigIIIsig, IIIpiIIIpi, IIIsIIIs, IIIsigIIIs, IIIsIIIsig;
/////////////////////////////////////

static const double f_sig=1.267;
static const double f_pi  =0.585; // these two factors are the orbital overlap fudge factors

///these are variables that i will use to make calculations faster/////////////////
static double EXP;
static double sin_h;
static double cos_h;
static double pt;

double p_pt[7]; // this will contain the powers of 1/pt^[i+1]
double p_p[7]; // this will contain the powers of 1/p^[i+1] -where i is the index!
///////////////////////////////////////////////////////////////////////////////////


static const double TOLL=1E-10; //for checking t ;0)

if (lbl1==0 && lbl2==0){
//	if(t<TOLL){
	   IsIs = _1S1S[P];
//	} // commented out this legacy code
/*	else{
		cout << "There is an error if you are here... this should only do H-H!" << endl;
		if(t<0) t=-t;
		pt = p*t;
		for(i=0;i<5;i++)p_pt[i]=pow(pt,-i-1);
		for(i=0;i<5;i++)p_p[i] =pow(p ,-i-1);
		cos_h = cosh(pt);
		sin_h = sinh(pt);
		EXP = exp(-p);
		IsIs = p*p*p/4*pow( (1-t*t), 1.5 ) * (A2(p_p, EXP) * B0(sin_h,p_pt) - B2(cos_h,sin_h,p_pt)*A0(p_p,EXP));
	//	out_pt_1s1s << pt << '\t'<< p << '\t' << IsIs <<endl;
	}*/

	F[c1][c2]=beta*IsIs;
}

else if (lbl1>1 && lbl2>1 &&lbl1<9 && lbl2<9 ){ //both row two
	if ( t < TOLL && t > -TOLL){ 
	        if ( lbl1 == 5 ){
			IIsIIs = _2S2S[P];

			IIsIIsig = _2S2Sig[P] ;

			IIsigIIsig = _2Sig2Sig[P];

			IIpiIIpi = _2Pi2Pi[P];
			IIsigIIs = IIsIIsig;

		}
		else {
		        EXP = exp(-p);
			p_p[0]=1.0/p;
			for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;

		    	IIsIIs     =  p*p*p*p*p/360.0 * (15.0 * A4( p_p, EXP) - 10.0 * A2( p_p, EXP) + 3.0 * A0( p_p, EXP) ) ;
			IIsIIsig   =  p*p*p*p*p / ( 60.0 * sqrt(3.0) ) * ( 5*A3(p_p, EXP) - A1( p_p, EXP) );
			IIsigIIsig =  p*p*p*p*p / 120.0 * ( 5.0 * A4(p_p, EXP) - 18.0 * A2(p_p, EXP) + 5.0 * A0(p_p, EXP) );
			IIpiIIpi   =  p*p*p*p*p / 120.0 * ( 5.0 * A4(p_p, EXP) -  6.0 * A2(p_p, EXP) + A0(p_p, EXP) );
			IIsigIIs   =  IIsIIsig;
		}
		//out_p_2s2s << p << '\t' << IIsIIs << endl;
		//out_p_2s2sig << p << '\t' <<  IIsIIsig <<endl;
	        //out_p_2sig2sig<< p<<'\t' << IIsigIIsig <<endl;
	        //out_p_2pi2pi << p <<'\t' << IIpiIIpi << endl;
	}
	else{

                pt = p*t;
//		for(i=0;i<5;i++)p_pt[i]=pow(pt,-i-1);
//		for(i=0;i<5;i++)p_p[i] =pow(p ,-i-1);

		p_p[0]=1.0/p;
		for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;

		p_pt[0]=1.0/pt;
		for(i=1;i<5;i++) p_pt[i] = p_pt[i-1] / pt;

		cos_h = cosh(pt);
		sin_h = sinh(pt);
		EXP = exp(-p);



		IIsIIs= p*p*p*p*p / 48.0 * pow( 1-t*t, 2.5) * ( A4(p_p,EXP) * B0(sin_h, p_pt) - 2.0 *
				A2(p_p, EXP) * B2(cos_h,sin_h,p_pt) + A0(p_p, EXP) * B4(cos_h, sin_h, p_pt) )  ;

                IIsIIsig= p*p*p*p*p / (16.0*sqrt(3.0)) * pow( 1-t*t, 2.5) * (
			  A3(p_p,EXP) * ( B0(sin_h, p_pt)      - B2(cos_h,sin_h,p_pt) ) +
			  A1(p_p,EXP) * ( B4(cos_h,sin_h,p_pt) - B2(cos_h,sin_h,p_pt) ) +
			  B1(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) - A4(p_p, EXP) ) 	+
			  B3(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) - A0(p_p, EXP) ) ) ;

                IIpiIIpi= p*p*p*p*p / 32.0 * pow( 1-t*t, 2.5) * (
			  A4(p_p, EXP) * ( B0(sin_h,p_pt) - B2(cos_h,sin_h,p_pt)) +
			  A2(p_p, EXP) * ( B4(cos_h,sin_h,p_pt) - B0(sin_h,p_pt)) +
			  A0(p_p, EXP) * ( B2(cos_h,sin_h,p_pt) - B4(cos_h,sin_h,p_pt)) );

                IIsigIIsig= p*p*p*p*p / 16.0 * pow( 1-t*t, 2.5) * (
			    B2(cos_h,sin_h,p_pt) * ( A0(p_p, EXP) + A4(p_p, EXP) ) -
			    A2(p_p, EXP) * ( B0(sin_h,p_pt) + B4(cos_h,sin_h,p_pt) ));
		for(i=0;i<5;i+=2) p_pt[i]=-p_pt[i];
		pt=-pt;
		cos_h = cosh(pt);
		sin_h = sinh(pt);
		IIsigIIs =  p*p*p*p*p / (16.0*sqrt(3.0)) * pow( 1-t*t, 2.5) * (
                          A3(p_p,EXP) * ( B0(sin_h, p_pt)      - B2(cos_h,sin_h,p_pt) ) +
                          A1(p_p,EXP) * ( B4(cos_h,sin_h,p_pt) - B2(cos_h,sin_h,p_pt) ) +
                          B1(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) - A4(p_p, EXP) )        +
                          B3(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) - A0(p_p, EXP) ) ) ;

		//out_pt_2s2s << pt << '\t'<< p <<'\t' << IIsIIs << endl;
		//out_pt_2s2sig << pt << '\t'<< p <<'\t' << IIsIIsig << endl;
		//out_pt_2sig2sig << pt << '\t'<< p <<'\t' << IIsigIIsig << endl;
		//out_pt_2pi2pi << pt << '\t'<< p <<'\t' << IIpiIIpi << endl;

	}

///////set the value of those elements that do not require calculation of geometric factors
        F[c1][c2]=beta * IIsIIs; // <2s|2s>
        F[c1+1][c2]=beta * IIsigIIs * (x); // <2px|2s>
        F[c1+2][c2]=beta * IIsigIIs * (y); // <2py|2s>
        F[c1+3][c2]=beta * IIsigIIs * (z); // <2pz|2s>

	F[c1][c2+1]=beta * IIsIIsig * (-x); // <2s| 2px>
	F[c1][c2+2]=beta * IIsIIsig * (-y);  //  <2s |2py> etc.
	F[c1][c2+3]=beta * IIsIIsig * (-z);
////////////////////////////////////////////////////////////////////////////////////////////

///////////in this section we calculate the geometric factors///////////
	        
 /*     general rotation about z by theta followed by one about x by psi

Rot[0][0] = cos(theta);
Rot[0][1] = cospsi * sin(theta);
Rot[0][2] = sin(psi) * sin(theta);   our x
Rot[1][0] = -sin(theta);
Rot[1][1] = cospsi * cos(theta);
Rot[1][2] = sin(psi) * cos(theta);   our y
Rot[2][1] = -sin(psi);
Rot[2][2] = cospsi;                our z

*/
        /* the rotation matrix that from unit z axis in the cartesian
         * coordinates to the axis x,y,z of the cylindrical coords
         * knowing the rotation matrix R, we know S in the cylindrical coord:
         *
         *             beta*f_pi*IIpiIIpi     0                 0
         *             0                     beta*f_pi*IIpiIIpi 0
         *             0                          0             beta*f_sig*IIsigIIsig
         *
         *
         *and we can determine S' by R S R'         
         
         */
        
        if ( z*z!= 1.) {
            double cospsi   = z;
            double sinpsi   = sqrt(1.-z*z);
            double costheta = y / sinpsi;
            double sintheta = x / sinpsi;
            double pi = IIpiIIpi * f_pi * beta;
            double sigma = IIsigIIsig * f_sig * beta;
            F[c1+1][c2+1] = costheta*costheta * pi + cospsi*cospsi * sintheta*sintheta * pi + sinpsi*sinpsi * sintheta*sintheta * sigma;
            F[c1+1][c2+2] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
            F[c1+1][c2+3] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
            F[c1+2][c2+1] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
            F[c1+2][c2+2] = sintheta*sintheta * pi + cospsi*cospsi * costheta*costheta * pi + sinpsi*sinpsi * costheta*costheta * sigma;
            F[c1+2][c2+3] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
            F[c1+3][c2+1] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
            F[c1+3][c2+2] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
            F[c1+3][c2+3] = sinpsi*sinpsi * pi + cospsi*cospsi * sigma;

            
        }
        else{
            double pi = IIpiIIpi * f_pi * beta;
            double sigma = IIsigIIsig * f_sig * beta;
            F[c1+1][c2+1] = sigma;
            F[c1+1][c2+2] = 0.;
            F[c1+1][c2+3] = 0.;
            F[c1+2][c2+1] = 0.;
            F[c1+2][c2+2] = pi;
            F[c1+2][c2+3] = 0.;
            F[c1+3][c2+1] = 0.;
            F[c1+3][c2+2] = 0.;
            F[c1+3][c2+3] = sigma;

        }	
        
}

else if (lbl1>1 && lbl1 <9 && lbl2==0){ //row 2 vs row 1
/*	if( t < TOLL &&  t > -TOLL){ // not going to happen really guv!
		EXP = exp(-p);
		p_p[0]=1.0/p;
		for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;

		IsIIs = p*p*p*p/( 12.0 * sqrt(3.0) ) * ( 3.0*A3(p_p,EXP) - A1(p_p,EXP) ) ; // could write sqrt(3)?
		IsIIsig = p*p*p*p/12.0*(3.0*A2(p_p,EXP) - A0(p_p,EXP) ) ;
		//out_p_1s2s << p << '\t' << IsIIs << endl;
		//out_p_1s2sig << p << '\t' << IsIIsig <<endl;
	}*/ 
//	else{ // commentin out these elses for speed
		t=-t; //Jo Dog! this line is WELL important, we can only calculate <1s|2s> , not <2sig|1s>: this latter overlap will be - the previous 
		if ( Abs(t + 0.15044247787610619469026548672566 ) < TOLL ){
		    IsIIs = _1S2S_hc[P];
		    IsIIsig = _1S2Sig_hc[P];
		}
		else {

			pt = p*t;
//		for(i=0;i<5;i++)p_pt[i]=pow(pt,-i-1);
//		for(i=0;i<5;i++)p_p[i] =pow(p ,-i-1);

			p_p[0]=1.0/p;
        	        for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;
			p_pt[0]=1.0/pt;
                	for(i=1;i<5;i++) p_pt[i] = p_pt[i-1] / pt;

			cos_h = cosh(pt);
			sin_h = sinh(pt);
			EXP = exp(-p);


			IsIIs =   p*p*p*p / ( 8.0 * sqrt(3.0 )) * pow(1-t, 2.5) * pow(1+t, 1.5) *
			( A3(p_p,EXP) * B0(sin_h, p_pt)        - A2(p_p, EXP) * B1(cos_h, sin_h, p_pt) -
			  A1(p_p,EXP) * B2(cos_h, sin_h, p_pt) + A0(p_p, EXP) * B3(cos_h, sin_h, p_pt) )  ;

			IsIIsig = p*p*p*p / 8.0 * pow(1-t, 2.5) * pow(1+t, 1.5) *
			(-A3(p_p,EXP) * B1(cos_h, sin_h, p_pt) + A2(p_p, EXP) * B0(sin_h, p_pt)       +
                          A1(p_p,EXP) * B3(cos_h, sin_h, p_pt) - A0(p_p, EXP) * B2(cos_h, sin_h, p_pt) ) ;
		}
//		cout << " A3: " << A3(p_p,EXP) << " B0 " << B0(sin_h, p_pt) << " A2 " << A2(p_p, EXP) << " B1: " << B1(cos_h, sin_h, p_pt) << " A1: " << A1(p_p,EXP) <<  " B2: " << B2(cos_h, sin_h, p_pt) << " A0: " << A0(p_p, EXP) << " B3: " << B3(cos_h, sin_h, p_pt) << endl;
       		//out_pt_1s2s  << pt << '\t'<< p <<'\t' << IsIIs << endl ;
		//out_pt_1s2sig << pt <<'\t'<< p <<'\t' << IsIIsig <<endl;
//	}

	F[c1][c2]=beta*IsIIs; // <2s|1s>
        F[c1+1][c2]=beta*(IsIIsig)*(x); // <2px|1s>
        F[c1+2][c2]=beta*(IsIIsig)*(y); // <2py|1s>
        F[c1+3][c2]=beta*(IsIIsig)*(z); // <2pz|1s>
}
else if (lbl1==0 && lbl2>1 && lbl2 <9){ //row 1 and row 2
//	cout << "nc2==1 && nc1==4"<< endl;
/*	if( t < TOLL &&  t > -TOLL){
                EXP = exp(-p);
                p_p[0]=1.0/p;
                for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;
                p_pt[0]=1.0/pt;
                for(i=1;i<5;i++) p_pt[i] = p_pt[i-1] / pt;


		IsIIs = p*p*p*p/( 12.0 * sqrt(3.0) ) * ( 3.0*A3(p_p,EXP) - A1(p_p,EXP) ); // could write sqrt(3)?
                IsIIsig = p*p*p*p/12.0*(3.0*A2(p_p,EXP) - A0(p_p,EXP) ) ;

		//out_p_1s2s << p << '\t' << IsIIs << endl;
		//out_p_1s2sig << p << '\t' << IsIIsig <<endl;
        }*/
//        else{
                if ( Abs(t + 0.1504424778761061946902654872566 ) < TOLL  ){
		   IsIIs = _1S2S_hc[P];
		   IsIIsig = _1S2Sig_hc[P];
		}
		else{
			pt = p*t;
//              for(i=0;i<5;i++)p_pt[i]=pow(pt,-i-1);
//              for(i=0;i<5;i++)p_p[i] =pow(p ,-i-1);
			p_pt[0]=1.0/pt;
			for(i=1;i<5;i++) p_pt[i] = p_pt[i-1] / pt;
			p_p[0]=1.0/p;
			for(i=1;i<5;i++) p_p[i] = p_p[i-1] / p;

        	        cos_h = cosh(pt);
                	sin_h = sinh(pt);
            		EXP = exp(-p);

              		IsIIs =   p*p*p*p / ( 8.0 * sqrt(3.0 )) * pow( 1.0 + t, 1.5) * pow (1.0 - t, 2.5) *
                        ( A3(p_p,EXP) * B0(sin_h, p_pt)        - A2(p_p, EXP) * B1(cos_h, sin_h, p_pt) -
                          A1(p_p,EXP) * B2(cos_h, sin_h, p_pt) + A0(p_p, EXP) * B3(cos_h, sin_h, p_pt) )  ;

                	IsIIsig = p*p*p*p / 8.0 * sqrt(1.0-t*t) * (1-t*t) * (1-t) *
                        (-A3(p_p,EXP) * B1(cos_h, sin_h, p_pt) + A2(p_p, EXP) * B0(sin_h, p_pt)       +
                          A1(p_p,EXP) * B3(cos_h, sin_h, p_pt) - A0(p_p, EXP) * B2(cos_h, sin_h, p_pt) ) ;
		}

//		cout << " A3: " << A3(p_p,EXP) << " B0 " << B0(sin_h, p_pt) << " A2 " << A2(p_p, EXP) << " B1: " << B1(cos_h, sin_h, p_pt) << " A1: " << A1(p_p,EXP) <<  " B2: " << B2(cos_h, sin_h, p_pt) << " A0: " << A0(p_p, EXP) << " B3: " << B3(cos_h, sin_h, p_pt) << endl;
		//out_pt_1s2s  << pt << '\t'<< p <<'\t' << IsIIs <<endl;
	        //out_pt_1s2sig << pt << '\t'<< p <<'\t' << IsIIsig <<endl;
  //      }

	F[c1][c2]=beta*IsIIs; //<1s|2s>
	F[c1][c2+1]=beta*IsIIsig*(-x);//<1s|2px>
	F[c1][c2+2]=beta*IsIIsig*(-y);//<1s|2py>
	F[c1][c2+3]=beta*IsIIsig*(-z);
}
else if(lbl1 > 9 && lbl2 > 9){ //row 3 row3
    if (t<TOLL && t > -TOLL){
	    EXP = exp(-p);	
	    p_p[0] =1.0/p;
	    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;
	    IIIsIIIs      = p*p*p*p*p*p*p/25200.0 *(35.0 * A6(p_p, EXP) -35.0 *A4(p_p,EXP) +21.0 * A2(p_p, EXP) -5.0 *A0(p_p,EXP) );
	    IIIsIIIsig    = p*p*p*p*p*p*p/(sqrt(3.0) * 12600.0) * (35.0*A5(p_p,EXP) -14.0*A3(p_p,EXP) +3.0*A1(p_p,EXP)   );
	    IIIsigIIIsig  = p*p*p*p*p*p*p/25200.0 *(35.0 * A6(p_p,EXP) -147.0 * A4(p_p, EXP) +85.0*A2(p_p, EXP)-21.0*A0(p_p,EXP) );
	    IIIpiIIIpi    = p*p*p*p*p*p*p/25200.0 *(35.0*A6(p_p,EXP) - 49.0*A4(p_p,EXP) + 17.0*A2(p_p, EXP) -3.0*A0(p_p,EXP));
	    IIIsigIIIs  = IIIsIIIsig;
    }
    else {
	    pt = p*t;
	    p_pt[0]=1.0/pt;
            for(i=1;i<7;i++) p_pt[i] = p_pt[i-1] / pt;
            cos_h = cosh(pt);
	    sin_h = sinh(pt);

	    EXP = exp(-p);	
	    p_p[0] =1.0/p;
	    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;

	    IIIsIIIs   =  p*p*p*p*p*p*p/1440.0 * (1-t)*(1-t)*(1-t)*sqrt(1-t) *
				(A6(p_p,EXP) *B0(cos_h,sin_h, p_pt) -3.0*A4(p_p,EXP) *B2(cos_h, sin_h, p_pt)
			       +3.0* A2(p_p,EXP) *B4(cos_h, sin_h, p_pt) - A0(p_p,EXP) *B6(cos_h,sin_h, p_pt)	)  ;
	    IIIsIIIsig = p*p*p*p*p*p*p/(sqrt(3.0) *480.0 ) * (1-t)*(1-t)*(1-t)*sqrt(1-t)*
	       		   (- A6(p_p, EXP)*B1(cos_h,sin_h,p_pt) + A5(p_p, EXP) *(B0(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt)) 
			    + A4(p_p, EXP)*(B1(cos_h,sin_h,p_pt)+2.0*B3(cos_h,sin_h,p_pt)) + 2.0 *A3(p_p, EXP) *(B4(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt))
			    - A2(p_p, EXP) * (2.0*B3(cos_h,sin_h,p_pt)+B5(cos_h,sin_h,p_pt)) + A1(p_p, EXP) * (B4(cos_h,sin_h,p_pt)-B6(cos_h,sin_h,p_pt)) 
			    + A0(p_p, EXP) * B5(cos_h,sin_h,p_pt));
	    IIIpiIIIpi   = p*p*p*p*p*p*p/(960.0 ) * (1-t)*(1-t)*(1-t)*sqrt(1-t) *
			   (A6(p_p, EXP) * (B0(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt)) + 
			    A4(p_p, EXP) * (2.0*B4(cos_h,sin_h,p_pt)-B0(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt)) +
			    A2(p_p, EXP) * (2.0*B2(cos_h,sin_h,p_pt)-B6(cos_h,sin_h,p_pt)-B4(cos_h,sin_h,p_pt)) +
			    A0(p_p, EXP) * (B6(cos_h,sin_h,p_pt)-B4(cos_h,sin_h,p_pt)) );
	    IIIsigIIIsig =  p*p*p*p*p*p*p/(480.0 ) * (1-t)*(1-t)*(1-t)*sqrt(1-t) * 
			    (A6(p_p, EXP) * B2(cos_h,sin_h,p_pt) - A4(p_p, EXP) *(B0(cos_h,sin_h,p_pt)+2.0*B4(cos_h,sin_h,p_pt)) 
			     + A2(p_p, EXP)*(B6(cos_h,sin_h,p_pt)+2.0*B2(cos_h,sin_h,p_pt)) - A0(p_p, EXP)*B4(cos_h,sin_h,p_pt) ) ;
	    t=-t;
            pt = p*t;
            p_pt[0]=1.0/pt;
            for(i=1;i<7;i+=2) p_pt[i] = -p_pt[i];
            cos_h = cosh(pt);
            sin_h = sinh(pt);
	    IIIsigIIIs   =  p*p*p*p*p*p*p/(sqrt(3.0) *480.0 ) * (1-t)*(1-t)*(1-t)*sqrt(1-t)*
	       		   (- A6(p_p, EXP)*B1(cos_h,sin_h,p_pt) + A5(p_p, EXP) *(B0(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt)) 
			    + A4(p_p, EXP)*(B1(cos_h,sin_h,p_pt)+2.0*B3(cos_h,sin_h,p_pt)) + 2.0 *A3(p_p, EXP) *(B4(cos_h,sin_h,p_pt)-B2(cos_h,sin_h,p_pt))
			    - A2(p_p, EXP) * (2.0*B3(cos_h,sin_h,p_pt)+B5(cos_h,sin_h,p_pt)) + A1(p_p, EXP) * (B4(cos_h,sin_h,p_pt)-B6(cos_h,sin_h,p_pt)) 
			    + A0(p_p, EXP) * B5(cos_h,sin_h,p_pt));
    }
    F[c1][c2]=beta * IIIsIIIs; // <3s|3s>
    F[c1+1][c2]=beta * IIIsIIIsig * (x); // <3px|3s>
    F[c1+2][c2]=beta * IIIsIIIsig * (y); // <3py|3s>
    F[c1+3][c2]=beta * IIIsIIIsig * (z); // <3pz|3s>

    F[c1][c2+1]=beta * IIIsigIIIs * (-x); // <3s| 3px>
    F[c1][c2+2]=beta * IIIsigIIIs * (-y);  //  <3s |3py> etc.
    F[c1][c2+3]=beta * IIIsigIIIs * (-z);
////////////////////////////////////////////////////////////////////////////////////////////

///////////in this section we calculate the geometric factors///////////
    if ( z*z!= 1.) {
        double cospsi   = z;
        double sinpsi   = sqrt(1.-z*z);
        double costheta = y / sinpsi;
        double sintheta = x / sinpsi;
        double pi = IIIpiIIIpi * f_pi * beta;
        double sigma = IIIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = costheta*costheta * pi + cospsi*cospsi * sintheta*sintheta * pi + sinpsi*sinpsi * sintheta*sintheta * sigma;
        F[c1+1][c2+2] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+1][c2+3] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+2][c2+1] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+2][c2+2] = sintheta*sintheta * pi + cospsi*cospsi * costheta*costheta * pi + sinpsi*sinpsi * costheta*costheta * sigma;
        F[c1+2][c2+3] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+1] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+3][c2+2] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+3] = sinpsi*sinpsi * pi + cospsi*cospsi * sigma;


    }
    else{
        double pi = IIIpiIIIpi * f_pi * beta;
        double sigma = IIIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = sigma;
        F[c1+1][c2+2] = 0.;
        F[c1+1][c2+3] = 0.;
        F[c1+2][c2+1] = 0.;
        F[c1+2][c2+2] = pi;
        F[c1+2][c2+3] = 0.;
        F[c1+3][c2+1] = 0.;
        F[c1+3][c2+2] = 0.;
        F[c1+3][c2+3] = sigma;

    }	
}

else if(lbl1 > 1 && lbl1 < 9 && lbl2 > 9){ //row 2 row 3
    //careful with Ts here!!!
    pt = p*t; 
    p_pt[0]=1.0/pt; 
    for(i=1;i<7;i++) p_pt[i] = p_pt[i-1] / pt; 
    cos_h = cosh(pt); 
    sin_h = sinh(pt); 
    EXP = exp(-p);	
    p_p[0] =1.0/p; 
    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;
    IIsIIIs     = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(30.0))*
       		  (A5(p_p, EXP) * B0(cos_h,sin_h,p_pt) - A4(p_p, EXP) *  B1(cos_h,sin_h,p_pt) - 2.0 * A3(p_p, EXP) * B2(cos_h,sin_h,p_pt) 
		   +2.0* A2(p_p, EXP) * B3(cos_h,sin_h,p_pt) + A1(p_p, EXP) * B4(cos_h,sin_h,p_pt) -A0(p_p, EXP) * B5(cos_h,sin_h,p_pt) );
    IIpiIIIpi   = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(32.0 * sqrt(30.0)) *
		  (A5(p_p, EXP) * (B0(cos_h,sin_h,p_pt) -B2(cos_h,sin_h,p_pt) ) + A4(p_p, EXP) * (B3(cos_h,sin_h,p_pt) - B1(cos_h,sin_h,p_pt)) + 
		   A3(p_p, EXP) * (B4(cos_h,sin_h,p_pt) -B0(cos_h,sin_h,p_pt))  + A2(p_p, EXP) * (B1(cos_h,sin_h,p_pt) - B5(cos_h,sin_h,p_pt)) + 
		   A1(p_p, EXP) * (B2(cos_h,sin_h,p_pt) -B4(cos_h,sin_h,p_pt) ) + A0(p_p, EXP) * (B5(cos_h,sin_h,p_pt) - B3(cos_h,sin_h,p_pt))  ) ;
    IIsIIIsig   = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(10.0))* 
		  (-A5(p_p, EXP) * B1(cos_h,sin_h,p_pt) + A4(p_p, EXP) * B0(cos_h,sin_h,p_pt) + 2.0 * A3(p_p, EXP) * B3(cos_h,sin_h,p_pt) - 
		   2.0 * A2(p_p, EXP) * B2(cos_h,sin_h,p_pt) - A1(p_p, EXP) * B5(cos_h,sin_h,p_pt) + A0(p_p, EXP) * B4(cos_h,sin_h,p_pt) ) ;
    IIsigIIIs  =  p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(10.0))* 
	          (A4(p_p, EXP) * (B0(cos_h,sin_h,p_pt) - 2.0 * B2(cos_h,sin_h,p_pt) ) + A1(p_p, EXP) * (2*B3(cos_h,sin_h,p_pt) - B5(cos_h,sin_h,p_pt) )  
		   + B1(cos_h,sin_h,p_pt) * (A5(p_p, EXP) - 2.0 * A3(p_p, EXP) ) +B4(cos_h,sin_h,p_pt) * (2.0*A2(p_p, EXP) - A0(p_p, EXP) ) );
    IIsigIIIsig = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(16.0*sqrt(30.0))*
       		  (A2(p_p, EXP) * (B1(cos_h,sin_h,p_pt) + B5(cos_h,sin_h,p_pt) ) - A3(p_p, EXP) * (B0(cos_h,sin_h,p_pt) + B4(cos_h,sin_h,p_pt) ) 
		   - B3(cos_h,sin_h,p_pt) * (A0(p_p, EXP) + A4(p_p, EXP) ) + B2(cos_h,sin_h,p_pt) * ( A1(p_p, EXP) +A5(p_p, EXP) )  ) ;
////////////////easy geometric factors
    F[c1][c2]=beta * IIsIIIs; // <2s|3s>
    F[c1+1][c2]=beta * IIsigIIIs * (x); // <2px|3s>
    F[c1+2][c2]=beta * IIsigIIIs * (y); // <2py|3s>
    F[c1+3][c2]=beta * IIsigIIIs * (z); // <2pz|3s>

    F[c1][c2+1]=beta * IIsIIIsig * (-x); // <2s| 3px>
    F[c1][c2+2]=beta * IIsIIIsig * (-y);  //  <2s |3py> etc.
    F[c1][c2+3]=beta * IIsIIIsig * (-z);
///////////in this section we calculate the geometric factors///////////
      if (z*z!= 1.) {
        double cospsi   = z;
        double sinpsi   = sqrt(1.-z*z);
        double costheta = y / sinpsi;
        double sintheta = x / sinpsi;
        double pi = IIpiIIIpi * f_pi * beta;
        double sigma = IIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = costheta*costheta * pi + cospsi*cospsi * sintheta*sintheta * pi + sinpsi*sinpsi * sintheta*sintheta * sigma;
        F[c1+1][c2+2] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+1][c2+3] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+2][c2+1] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+2][c2+2] = sintheta*sintheta * pi + cospsi*cospsi * costheta*costheta * pi + sinpsi*sinpsi * costheta*costheta * sigma;
        F[c1+2][c2+3] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+1] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+3][c2+2] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+3] = sinpsi*sinpsi * pi + cospsi*cospsi * sigma;


    }
    else{
        double pi = IIpiIIIpi * f_pi * beta;
        double sigma = IIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = sigma;
        F[c1+1][c2+2] = 0.;
        F[c1+1][c2+3] = 0.;
        F[c1+2][c2+1] = 0.;
        F[c1+2][c2+2] = pi;
        F[c1+2][c2+3] = 0.;
        F[c1+3][c2+1] = 0.;
        F[c1+3][c2+2] = 0.;
        F[c1+3][c2+3] = sigma;

    }	
}
else if( lbl1 > 9 && lbl2 > 1 && lbl2 < 9 ){ //row 3 row 2
//careful with Ts here!!!
   //  cout << "Warning! Using untested Sulphur code! Take care!" <<endl;   
    t = -t;//check this out!
    pt = p*t; 
    p_pt[0]=1.0/pt; 
    for(i=1;i<7;i++) p_pt[i] = p_pt[i-1] / pt; 
    cos_h = cosh(pt); 
    sin_h = sinh(pt); 
    EXP = exp(-p);	
    p_p[0] =1.0/p; 
    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;
    IIsIIIs     = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(30.0))*
       		  (A5(p_p, EXP) * B0(cos_h,sin_h,p_pt) - A4(p_p, EXP) *  B1(cos_h,sin_h,p_pt) - 2.0 * A3(p_p, EXP) * B2(cos_h,sin_h,p_pt) 
		   +2.0*A2(p_p, EXP) * B3(cos_h,sin_h,p_pt) + A1(p_p, EXP) * B4(cos_h,sin_h,p_pt) -A0(p_p, EXP) * B5(cos_h,sin_h,p_pt) );
    IIpiIIIpi   = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(32.0 * sqrt(30.0)) *
		  (A5(p_p, EXP) * (B0(cos_h,sin_h,p_pt) -B2(cos_h,sin_h,p_pt) ) + A4(p_p, EXP) * (B3(cos_h,sin_h,p_pt) - B1(cos_h,sin_h,p_pt)) + 
		   A3(p_p, EXP) * (B4(cos_h,sin_h,p_pt) -B0(cos_h,sin_h,p_pt))  + A2(p_p, EXP) * (B1(cos_h,sin_h,p_pt) - B5(cos_h,sin_h,p_pt)) + 
		   A1(p_p, EXP) * (B2(cos_h,sin_h,p_pt) -B4(cos_h,sin_h,p_pt) ) + A0(p_p, EXP) * (B5(cos_h,sin_h,p_pt) - B3(cos_h,sin_h,p_pt))  ) ;
    IIIsigIIs   = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(10.0))* 
		  (-A5(p_p, EXP) * B1(cos_h,sin_h,p_pt) + A4(p_p, EXP) * B0(cos_h,sin_h,p_pt) + 2.0 * A3(p_p, EXP) * B3(cos_h,sin_h,p_pt) - 
		   2.0 * A2(p_p, EXP) * B2(cos_h,sin_h,p_pt) - A1(p_p, EXP) * B5(cos_h,sin_h,p_pt) + A0(p_p, EXP) * B4(cos_h,sin_h,p_pt) ) ;
    IIIsIIsig  =  p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(48.0*sqrt(10.0))* 
	          (A4(p_p, EXP) * (B0(cos_h,sin_h,p_pt) - 2.0 * B2(cos_h,sin_h,p_pt) ) + A1(p_p, EXP) * (2*B3(cos_h,sin_h,p_pt) - B5(cos_h,sin_h,p_pt) )  
		   + B1(cos_h,sin_h,p_pt) * (A5(p_p, EXP) - 2.0 * A3(p_p, EXP) ) +B4(cos_h,sin_h,p_pt) * (2.0*A2(p_p, EXP) - A0(p_p, EXP) ) );
    IIsigIIIsig = p*p*p*p*p*p*pow(1+t,2.5)*pow(1-t,3.5)/(16.0*sqrt(30.0))*
       		  (A2(p_p, EXP) * (B1(cos_h,sin_h,p_pt) + B5(cos_h,sin_h,p_pt) ) - A3(p_p, EXP) * (B0(cos_h,sin_h,p_pt) + B4(cos_h,sin_h,p_pt) ) 
		   - B3(cos_h,sin_h,p_pt) * (A0(p_p, EXP) + A4(p_p, EXP) ) + B2(cos_h,sin_h,p_pt) * ( A1(p_p, EXP) +A5(p_p, EXP) )  ) ;
////////////////easy geometric factors (might be the other wat around)
    F[c1][c2]=beta * IIsIIIs; // <2s|3s>
    F[c1+1][c2]=beta * IIIsigIIs * (x); // <3px|2s>
    F[c1+2][c2]=beta * IIIsigIIs * (y); // <3py|2s>
    F[c1+3][c2]=beta * IIIsigIIs * (z); // <3pz|2s>

    F[c1][c2+1]=beta * IIIsIIsig * (-x); // <3s| 2px>
    F[c1][c2+2]=beta * IIIsIIsig * (-y);  //  <3s |2py> etc.
    F[c1][c2+3]=beta * IIIsIIsig * (-z);
    
///////////in this section we calculate the geometric factors///////////
    if ( z*z!= 1.) {
        double cospsi   = z;
        double sinpsi   = sqrt(1.-z*z);
        double costheta = y / sinpsi;
        double sintheta = x / sinpsi;
        double pi = IIpiIIIpi * f_pi * beta;
        double sigma = IIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = costheta*costheta * pi + cospsi*cospsi * sintheta*sintheta * pi + sinpsi*sinpsi * sintheta*sintheta * sigma;
        F[c1+1][c2+2] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+1][c2+3] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+2][c2+1] = -costheta * sintheta * pi + cospsi*cospsi * sintheta * costheta * pi + sinpsi*sinpsi * sintheta * costheta * sigma;
        F[c1+2][c2+2] = sintheta*sintheta * pi + cospsi*cospsi * costheta*costheta * pi + sinpsi*sinpsi * costheta*costheta * sigma;
        F[c1+2][c2+3] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+1] = -cospsi * sintheta * sinpsi * pi + sinpsi * sintheta * cospsi * sigma;
        F[c1+3][c2+2] = -cospsi * costheta * sinpsi * pi + sinpsi * costheta * cospsi * sigma;
        F[c1+3][c2+3] = sinpsi*sinpsi * pi + cospsi*cospsi * sigma;


    }
    else{
        double pi = IIpiIIIpi * f_pi * beta;
        double sigma = IIsigIIIsig * f_sig * beta;
        F[c1+1][c2+1] = sigma;
        F[c1+1][c2+2] = 0.;
        F[c1+1][c2+3] = 0.;
        F[c1+2][c2+1] = 0.;
        F[c1+2][c2+2] = pi;
        F[c1+2][c2+3] = 0.;
        F[c1+3][c2+1] = 0.;
        F[c1+3][c2+2] = 0.;
        F[c1+3][c2+3] = sigma;

    }	
}
else if( lbl1 <2 && lbl2 >9){ //row 1 row 3
    pt = p*t; 
    p_pt[0]=1.0/pt; 
    for(i=1;i<7;i++) p_pt[i] = p_pt[i-1] / pt; 
    cos_h = cosh(pt); 
    sin_h = sinh(pt); 
    EXP = exp(-p);	
    p_p[0] =1.0/p; 
    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;

    IsIIIs   = p*p*p*p*p*pow(1+t,1.5 )*pow(1-t,3.5 )/(24.0 * sqrt(10.0) )*
       		( A4(p_p, EXP) * B0(cos_h,sin_h,p_pt) -2.0 * A3(p_p, EXP) * B1(cos_h,sin_h,p_pt) 
		  + 2.0 * A1(p_p, EXP) * B3(cos_h,sin_h,p_pt) -A0(p_p, EXP) * B4(cos_h,sin_h,p_pt) )	;
    IsIIIsig = p*p*p*p*p*pow(1+t,1.5 )*pow(1-t,3.5 )/(8.0 * sqrt(30.0) )*
		( A3(p_p, EXP) * (B0(cos_h,sin_h,p_pt) + B2(cos_h,sin_h,p_pt) ) - A1(p_p, EXP) * ( B2(cos_h,sin_h,p_pt) + B4(cos_h,sin_h,p_pt) )
		  - B1(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) + A4(p_p, EXP)) + B3(cos_h,sin_h,p_pt) * ( A0(p_p, EXP) + A2(p_p, EXP) ) );
 
    F[c1][c2]=beta*IsIIs; //<1s|3s>
    F[c1][c2+1]=beta*IsIIsig*(-x);//<1s|3px>
    F[c1][c2+2]=beta*IsIIsig*(-y);//<1s|3py>
    F[c1][c2+3]=beta*IsIIsig*(-z);
}
else if(lbl1 > 9 && lbl2 < 1){ //row 3 row 1
    t=-t;
    pt = p*t; 
    p_pt[0]=1.0/pt; 
    for(i=1;i<7;i++) p_pt[i] = p_pt[i-1] / pt; 
    cos_h = cosh(pt); 
    sin_h = sinh(pt); 
    EXP = exp(-p);	
    p_p[0] =1.0/p; 
    for (i =1;i<7;i++) p_p[i] = p_p[i-1]/p;
 
    IsIIIs   = p*p*p*p*p*pow(1+t,1.5 )*pow(1-t,3.5 )/(24.0 * sqrt(10.0) )*
       		( A4(p_p, EXP) * B0(cos_h,sin_h,p_pt) -2.0 * A3(p_p, EXP) * B1(cos_h,sin_h,p_pt) 
		  + 2.0 * A1(p_p, EXP) * B3(cos_h,sin_h,p_pt) -A0(p_p, EXP) * B4(cos_h,sin_h,p_pt) )	;
    IsIIIsig = p*p*p*p*p*pow(1+t,1.5 )*pow(1-t,3.5 )/(8.0 * sqrt(30.0) )*
		( A3(p_p, EXP) * (B0(cos_h,sin_h,p_pt) + B2(cos_h,sin_h,p_pt) ) - A1(p_p, EXP) * ( B2(cos_h,sin_h,p_pt) + B4(cos_h,sin_h,p_pt) )
		  - B1(cos_h,sin_h,p_pt) * ( A2(p_p, EXP) + A4(p_p, EXP)) + B3(cos_h,sin_h,p_pt) * ( A0(p_p, EXP) + A2(p_p, EXP) ) );
 
    F[c1][c2]=beta*IsIIIs; // <3s|1s>
    F[c1+1][c2]=beta*(IsIIIsig)*(x); // <3px|1s>
    F[c1+2][c2]=beta*(IsIIIsig)*(y); // <3py|1s>
    F[c1+3][c2]=beta*(IsIIIsig)*(z); // <3pz|1s>
}
else{
	cerr << "error in the input of find F function, expecting only row 1 , row 2 or row 3 elements" << endl;
}

}

void fock::calc_F_lean( ) const { 

int n_mol1 = molecules.first  -> getN();
int n_mol2 = molecules.second -> getN();

int i,j;

///////////////////////////////////////////////
////////generate Fi////////////////////////////


// fill up F, one element - i.e. more than 1 basis set at a time
int nb1_tmp; // counts how many basis sets we need to do at a time
int nb2_tmp;

double beta_av; // the value of the average beta
double p; // these two variables have the same meaning as in the mulliken paper
double t;

int b_c1; // counter for the basis set on mol1 that we are doing
int b_c2; // counter for the basis set on mol2 that we are doing

int a,b; // temporary counters to make writing easier

double dx,dy,dz; // doubles describing the relative position of 2 and 1
double length;

b_c1=0;
for(i=0;i<n_mol1;i++){

	a=  (molecules.first) -> getlbl(i);
        if (a!=18){ //KIND OF BOTCHY... skipping the ghost atoms basically

            if(a<=2) nb1_tmp=1;
            else nb1_tmp=4; 
            b_c2=0; // set the counter for mol2 to zero
            for(j=0;j<n_mol2;j++) {

                    b= (molecules.second) -> getlbl(j);
                    if (b!=18){
        ///////////set common position///////////////////////
                        dx = (molecules.second)->getx(j) -  (molecules.first)->getx(i) ;
                        dy = (molecules.second)->gety(j) -  (molecules.first)->gety(i) ;
                        dz = (molecules.second)->getz(j) -  (molecules.first)->getz(i) ;

                        length = sqrt( dx*dx + dy*dy + dz*dz);
                        dx = dx / length;
                        dy = dy / length;
                        dz = dz / length;
        //		length = length / Ra; // express the length in Bohr
        //////////////////////////////////////////////////////

        ///////////set t and p///////////////////////////////

                        p = ( Mu[a] + Mu[b] ) * length / 2;
                        t = ( Mu[a] - Mu[b] ) / ( Mu[a] + Mu[b] ) ; // in the forMulae on the paper t => 0
        /////////////////////////////////////////////////////

        ///////////calculate av_beta//////////////////////////////
                        beta_av = (Beta[a] + Beta[b])/2.0;
        //////////////////////////////////////////////////////////

                        if (b<2) nb2_tmp=1;
                        else nb2_tmp=4;


                        calc_F_el_with_HASH(p,t,beta_av,dx,dy,dz,b_c1,b_c2,a,b);
                        b_c2+=nb2_tmp;
                    }
            }
            b_c1+=nb1_tmp;
        }
}
}

vector <double> fock::calcJ ( vector <pair <unsigned int, unsigned int> > input) const {
	calc_F_lean();
	vector < pair<unsigned int, unsigned int> > :: iterator iter;
	vector <double> res;
	int nbasis_1 = molecules.first -> getNBasis();
	int nbasis_2 = molecules.second-> getNBasis();
	for ( iter = input.begin() ; iter < input.end() ; ++ iter ) {
		double *orb1 = molecules.first -> getorb(iter->first);
		double *orb2 = molecules.second-> getorb(iter->second);	
		double J =0.;
		double *p_F = &F[0][0];
		for ( int i = 0; i <nbasis_1 ; ++i ) {
			for ( int j =0 ; j<nbasis_2 ;++ j ){
			 	J += (*p_F) * orb1[i] * orb2[j];
				++p_F;
			}
		}
		res.push_back(J);
	}
	return res;
}

double fock::calcJ ( pair <int, int> input) const {
	calc_F_lean();

	int nbasis_1 = molecules.first -> getNBasis();
	int nbasis_2 = molecules.second-> getNBasis();
	double *orb1 = molecules.first -> getorb(input.first);
	double *orb2 = molecules.second-> getorb(input.second);	
	double J =0.;
	double *p_F = &F[0][0];
	for ( int i = 0; i <nbasis_1 ; ++i ) {
		for ( int j =0 ; j<nbasis_2 ;++ j ){
		    J += (*p_F) * orb1[i] * orb2[j];
		    ++p_F;
	    }
	}
	
	return J;
}

}}
