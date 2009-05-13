#include "orbitals.h"
#include <stdio.h>
#include <stdlib.h>


int (orb::*orb::read_orb)(const char *)=&orb::read_orb_gauss;

void orb::init_orbitals (string * basis, const int & N, const char * namefile ){
    NBasis = N;
    bs = new string [NBasis];
    psi = new double* [NBasis];
    psi[0] = new double [NBasis * NBasis];
    bs[0]=basis[0];
    for ( int i = 1 ; i < NBasis ; i ++){
            bs[i] = basis[i];
            psi[i] = psi[i-1] + NBasis;
    }
    (this->*read_orb) ( namefile);
    /*reorder the d orbitals to agree with libint */
    int count=0;
    while (count < NBasis ){
    	if (bs[count] == "s") count +=1;
	else if (bs[count] == "x") count +=3;
	else if (bs[count] == "xx") { 
		bs[count+1] = "xy";
		bs[count+2] = "xz";
		bs[count+3] = "yy";	
		bs[count+4] = "yz";
		bs[count+5] = "zz";
		for (int j=0; j<NBasis; ++j){
		    double tyy = psi[j][count+1];
		    double tzz = psi[j][count+2];
		    double txy = psi[j][count+3];
		    double txz = psi[j][count+4];
		    double tyz = psi[j][count+5];

		    psi[j][count+1] = txy; 
		    psi[j][count+2] = txz;
		    psi[j][count+3] = tyy;
		    psi[j][count+4] = tyz; 
		    psi[j][count+5] = tzz;
		}
		count +=6;
	}
	else {
	    cerr << bs[count] <<  " ERROR in reordering the orbitals" <<endl;
	}
    }
}

void orb::strip_orbitals (const vector <unsigned int>& a){
    int nrorbs = a.size();
    double** psinew = new double*[nrorbs];
    psinew[0] = new double [nrorbs*NBasis];
    for(int i=1; i<nrorbs; i++){
        psinew[i] = psinew[i-1]+NBasis;
    }
    for (int i=0; i<nrorbs; i++){
        for (int j=0; j<NBasis; j++){
            psinew[i][j] = psi[a[i]][j];
        }
    }
    delete [] psi[0];
    delete [] psi;
    psi = new double *[nrorbs];
    psi[0] = new double [nrorbs*NBasis];
    for(int i=1; i<nrorbs; i++){
        psi[i] = psi[i-1]+NBasis;
    }
    for(int i=0; i<nrorbs; i++){
        for (int j=0; j<NBasis; j++){
            psi[i][j]=psinew[i][j];
        }
    }
    #ifdef DEBUG
    cout << "nrorbs: " << nrorbs << endl;
    cout << "psi[0][i] :" << endl;
    for (int j=0; j<NBasis; j++){
        cout << psi[0][j] << endl;
    }
    #endif
}

void orb::init_orbitals_stripped(const orb& orb1, const int& nrorbs ){
    NBasis = orb1.NBasis;
    bs = new string [NBasis];
    psi = new double* [nrorbs];
    psi[0] = new double [nrorbs*NBasis];
    bs[0]=orb1.bs[0];
    for(int i=1; i<nrorbs; i++){
        psi[i] = psi[i-1]+NBasis;
    }
    for(int i = 0; i < nrorbs; i++){
        for ( int j = 0 ; j < NBasis ; j ++){
            bs[i] = orb1.bs[i];
            psi[i][j] = orb1.psi[i][j];
            //cout << "(" << i << "," << j << "): " << psi[i][j] << endl;
        }
    }
}

void orb::set_read_orb_gam(){
    read_orb=&orb::read_orb_gamess;
}

void orb::set_read_orb_gauss(){
    read_orb=&orb::read_orb_gauss;
}

void orb::rot_orb(const int & orb , const double rot[3][3] ){
	int i=0;
	while (i<NBasis)
	{
		if ( bs[i] == "s"){
			i++;
		}

		else if (bs[i] == "x" ){	
			double x_t,y_t,z_t;
			x_t = psi[orb][i]*rot[0][0] + psi[orb][i+1] * rot[0][1] + psi[orb][i+2] * rot[0][2];
			y_t = psi[orb][i]*rot[1][0] + psi[orb][i+1] * rot[1][1] + psi[orb][i+2] * rot[1][2];
			z_t = psi[orb][i]*rot[2][0] + psi[orb][i+1] * rot[2][1] + psi[orb][i+2] * rot[2][2] ;
			psi[orb][i]=x_t;
			psi[orb][i+1]=y_t;
			psi[orb][i+2]=z_t;
			i=i+3;
		}
		else if (bs[i] =="xx" ){
		    	double xx_t, xy_t, xz_t, yy_t, yz_t, zz_t;
			xx_t= rot[0][0]*rot[0][0] * psi[orb][i]   + rot[0][0]*rot[0][1] * psi[orb][i+1] + rot[0][0]*rot[0][2]*psi[orb][i+2]
			    + rot[0][1]*rot[0][0] * psi[orb][i+1] + rot[0][1]*rot[0][1] * psi[orb][i+3] + rot[0][1]*rot[0][2]*psi[orb][i+4]
			    + rot[0][2]*rot[0][0] * psi[orb][i+2] + rot[0][2]*rot[0][1] * psi[orb][i+4] + rot[0][2]*rot[0][2]*psi[orb][i+5]
			    ;
			xy_t= rot[0][0]*rot[1][0] * psi[orb][i]   + rot[0][0]*rot[1][1] * psi[orb][i+1] + rot[0][0]*rot[1][2]*psi[orb][i+2]
			    + rot[0][1]*rot[1][0] * psi[orb][i+1] + rot[0][1]*rot[1][1] * psi[orb][i+3] + rot[0][1]*rot[1][2]*psi[orb][i+4]
			    + rot[0][2]*rot[1][0] * psi[orb][i+2] + rot[0][2]*rot[1][1] * psi[orb][i+4] + rot[0][2]*rot[1][2]*psi[orb][i+5]
			    ;
			xz_t= rot[0][0]*rot[2][0] * psi[orb][i]   + rot[0][0]*rot[2][1] * psi[orb][i+1] + rot[0][0]*rot[2][2]*psi[orb][i+2]
			    + rot[0][1]*rot[2][0] * psi[orb][i+1] + rot[0][1]*rot[2][1] * psi[orb][i+3] + rot[0][1]*rot[2][2]*psi[orb][i+4]
			    + rot[0][2]*rot[2][0] * psi[orb][i+2] + rot[0][2]*rot[2][1] * psi[orb][i+4] + rot[0][2]*rot[2][2]*psi[orb][i+5]
			    ;
			yy_t= rot[1][0]*rot[1][0] * psi[orb][i]   + rot[1][0]*rot[1][1] * psi[orb][i+1] + rot[1][0]*rot[1][2]*psi[orb][i+2]
			    + rot[1][1]*rot[1][0] * psi[orb][i+1] + rot[1][1]*rot[1][1] * psi[orb][i+3] + rot[1][1]*rot[1][2]*psi[orb][i+4]
			    + rot[1][2]*rot[1][0] * psi[orb][i+2] + rot[1][2]*rot[1][1] * psi[orb][i+4] + rot[1][2]*rot[1][2]*psi[orb][i+5]
			    ;
			yz_t= rot[1][0]*rot[2][0] * psi[orb][i]   + rot[1][0]*rot[2][1] * psi[orb][i+1] + rot[1][0]*rot[2][2]*psi[orb][i+2]
			    + rot[1][1]*rot[2][0] * psi[orb][i+1] + rot[1][1]*rot[2][1] * psi[orb][i+3] + rot[1][1]*rot[2][2]*psi[orb][i+4]
			    + rot[1][2]*rot[2][0] * psi[orb][i+2] + rot[1][2]*rot[2][1] * psi[orb][i+4] + rot[1][2]*rot[2][2]*psi[orb][i+5]
			    ;
			zz_t= rot[2][0]*rot[2][0] * psi[orb][i]   + rot[2][0]*rot[2][1] * psi[orb][i+1] + rot[2][0]*rot[2][2]*psi[orb][i+2]
			    + rot[2][1]*rot[2][0] * psi[orb][i+1] + rot[2][1]*rot[2][1] * psi[orb][i+3] + rot[2][1]*rot[2][2]*psi[orb][i+4]
			    + rot[2][2]*rot[2][0] * psi[orb][i+2] + rot[2][2]*rot[2][1] * psi[orb][i+4] + rot[2][2]*rot[2][2]*psi[orb][i+5]
			    ;

			psi[orb][i]   = xx_t;
			psi[orb][i+1] = xy_t;
			psi[orb][i+2] = xz_t;
			psi[orb][i+3] = yy_t;
			psi[orb][i+4] = yz_t;
			psi[orb][i+5] = zz_t;

			i+=6;
			    
		}
		else {
			cerr << " Error in  rot about axis" <<endl;
		}
	}
}

inline void orb::rot_orb(const double rot[3][3]){
	for (int i =0; i< NBasis;i++){
		rot_orb( i, rot);
	}
}

/*void orb::rot_orbs(const vector <unsigned int>& orbs, int* i, double* psi2, const matrix& rot){
    int maxcount = orbs.size();
    for (int count=0; count<maxcount; count++){
        rot_orb(count, i, psi2, rot);
    }
}*/

void orb::rot_orb(const int & orb , int *i, double * psi2, const matrix &rot){
    //cout << "BS[" << *i << "] = " << bs[*i] << endl;
    //cout << "psi[orb][i] = psi[" << orb << "][" << *i <<  "]:" << psi[orb][*i] << endl;
    
    
            if ( bs[*i] == "s"){
                    (*i)++;
            }
            else if (bs[*i] == "x" ){	
                    double x_t,y_t,z_t;
                    psi2[*i] = psi[orb][*i]*rot.get(0,0) + psi[orb][*i+1] * rot.get(0,1) + psi[orb][*i+2] * rot.get(0,2);
                    psi2[*i+1]= psi[orb][*i]*rot.get(1,0) + psi[orb][*i+1] * rot.get(1,1) + psi[orb][*i+2] * rot.get(1,2);
                    psi2[*i+2]= psi[orb][*i]*rot.get(2,0) + psi[orb][*i+1] * rot.get(2,1) + psi[orb][*i+2] * rot.get(2,2) ;
                    (*i)+=3;
            }
            else if (bs[*i] =="xx" ){
                    double xx_t, xy_t, xz_t, yy_t, yz_t, zz_t;
                    psi2[*i]= rot.get(0,0)*rot.get(0,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(0,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(0,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(0,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(0,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(0,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(0,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(0,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(0,2)*psi[orb][*i+5]
                        ;
                    psi2[*i+1]= rot.get(0,0)*rot.get(1,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(1,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(1,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(1,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(1,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(1,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(1,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(1,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(1,2)*psi[orb][*i+5]
                        ;
                    psi2[*i+2]= rot.get(0,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;
                    psi2[*i+3]= rot.get(1,0)*rot.get(1,0) * psi[orb][*i]   + rot.get(1,0)*rot.get(1,1) * psi[orb][*i+1] + rot.get(1,0)*rot.get(1,2)*psi[orb][*i+2]
                        + rot.get(1,1)*rot.get(1,0) * psi[orb][*i+1] + rot.get(1,1)*rot.get(1,1) * psi[orb][*i+3] + rot.get(1,1)*rot.get(1,2)*psi[orb][*i+4]
                        + rot.get(1,2)*rot.get(1,0) * psi[orb][*i+2] + rot.get(1,2)*rot.get(1,1) * psi[orb][*i+4] + rot.get(1,2)*rot.get(1,2)*psi[orb][*i+5]
                        ;
                    psi2[*i+4]= rot.get(1,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(1,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(1,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(1,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(1,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(1,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(1,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(1,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(1,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;
                    psi2[*i+5]= rot.get(2,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(2,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(2,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(2,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(2,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(2,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(2,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(2,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(2,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;

                    (*i)+=6;

            }
            else {
                    cerr << " Error in rot about axis" <<endl;
            }
}

void orb::rot_orb(const int & orb , int *i, const matrix &rot){

            if ( bs[*i] == "s"){
                    (*i)++;
            }

            else if (bs[*i] == "x" ){	
                    double x_t,y_t,z_t;
                    x_t = psi[orb][*i]*rot.get(0,0) + psi[orb][*i+1] * rot.get(0,1) + psi[orb][*i+2] * rot.get(0,2);
                    y_t = psi[orb][*i]*rot.get(1,0) + psi[orb][*i+1] * rot.get(1,1) + psi[orb][*i+2] * rot.get(1,2);
                    z_t = psi[orb][*i]*rot.get(2,0) + psi[orb][*i+1] * rot.get(2,1) + psi[orb][*i+2] * rot.get(2,2) ;
                    psi[orb][*i]=x_t;
                    psi[orb][*i+1]=y_t;
                    psi[orb][*i+2]=z_t;
                    (*i)+=3;
            }
            else if (bs[*i] =="xx" ){
                    double xx_t, xy_t, xz_t, yy_t, yz_t, zz_t;
                    xx_t= rot.get(0,0)*rot.get(0,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(0,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(0,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(0,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(0,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(0,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(0,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(0,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(0,2)*psi[orb][*i+5]
                        ;
                    xy_t= rot.get(0,0)*rot.get(1,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(1,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(1,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(1,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(1,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(1,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(1,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(1,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(1,2)*psi[orb][*i+5]
                        ;
                    xz_t= rot.get(0,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(0,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(0,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(0,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(0,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(0,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(0,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(0,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(0,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;
                    yy_t= rot.get(1,0)*rot.get(1,0) * psi[orb][*i]   + rot.get(1,0)*rot.get(1,1) * psi[orb][*i+1] + rot.get(1,0)*rot.get(1,2)*psi[orb][*i+2]
                        + rot.get(1,1)*rot.get(1,0) * psi[orb][*i+1] + rot.get(1,1)*rot.get(1,1) * psi[orb][*i+3] + rot.get(1,1)*rot.get(1,2)*psi[orb][*i+4]
                        + rot.get(1,2)*rot.get(1,0) * psi[orb][*i+2] + rot.get(1,2)*rot.get(1,1) * psi[orb][*i+4] + rot.get(1,2)*rot.get(1,2)*psi[orb][*i+5]
                        ;
                    yz_t= rot.get(1,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(1,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(1,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(1,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(1,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(1,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(1,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(1,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(1,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;
                    zz_t= rot.get(2,0)*rot.get(2,0) * psi[orb][*i]   + rot.get(2,0)*rot.get(2,1) * psi[orb][*i+1] + rot.get(2,0)*rot.get(2,2)*psi[orb][*i+2]
                        + rot.get(2,1)*rot.get(2,0) * psi[orb][*i+1] + rot.get(2,1)*rot.get(2,1) * psi[orb][*i+3] + rot.get(2,1)*rot.get(2,2)*psi[orb][*i+4]
                        + rot.get(2,2)*rot.get(2,0) * psi[orb][*i+2] + rot.get(2,2)*rot.get(2,1) * psi[orb][*i+4] + rot.get(2,2)*rot.get(2,2)*psi[orb][*i+5]
                        ;

                    psi[orb][*i]   = xx_t;
                    psi[orb][*i+1] = xy_t;
                    psi[orb][*i+2] = xz_t;
                    psi[orb][*i+3] = yy_t;
                    psi[orb][*i+4] = yz_t;
                    psi[orb][*i+5] = zz_t;

                    (*i)+=6;

            }
            else {
                    cerr << " Error in  rot about axis" <<endl;
            }
}


void orb::rotate_someatoms(vector<int> atoms , matrix * M, 
        double * psi2 , const int & _orb){
        //cout << "orb = " << _orb << endl;
        vector <int>::iterator it;
        for (it = atoms.begin() ; it != atoms.end(); ++it){
            int first_basis = (_basis_on_atom[*it]).first;
            int last_basis = first_basis + (_basis_on_atom[*it]).second;
            int i=first_basis;
            while (i < last_basis){
                rot_orb(_orb, &i, psi2, *M); /// JAMES HERE FOR FUCKS SAKE!!!
            }
        }
}

/*void orb::rotate_someatoms(vector <int> atoms, matrix* M, double *psi2, const vector <unsigned int>& _orbs){
   vector <int>::iterator it;
   for (it = atoms.begin(); it != atoms.end(); ++it){
        int first_basis = (_basis_on_atom[*it]).first;
        int last_basis = first_basis + (_basis_on_atom[*it]).second;
        int i=first_basis;
        while (i < last_basis){
            rot_orbs(_orbs, &i, psi2, *M);
        }
    }
}*/

void orb::rot_orb(const int & orb , const matrix &rot){
	int i=0;
	while (i<NBasis)
	{
		rot_orb(orb, &i , rot);
	}
}

inline void orb::rot_orb(const matrix &rot){
	for (int i =0; i< NBasis;i++){
		rot_orb( i, rot);
	}
}


int orb::read_orb_gauss( const char * nameorbs)
{
 ifstream in(nameorbs); // read in from the average file
 if (!in){
     cout << "Error file  nameorbs  containg the orbitals non existant " <<endl;
     return 1;
 }

 string line;
 string number;

 int i,j;

 vector <string> file;

 int n_lin_in_wt = (NBasis-1)/5+2;
 int k=0;

 //skip the fortran specification line
 getline (in, line);
 while (in){
	getline (in, line);
        if(k  % (n_lin_in_wt) != 0)
        {
                while (line.size() > 1 )
                {
                        number.assign(line,0,15);
                        file.push_back(number);
                        line.erase(0,15);
                }
        }
        k++;
   }

   for(i=0,j=0,k=0;i<file.size();i++,k++)
   {

        file[i].replace(file[i].size()-4,1,1,'e');
        sscanf(file[i].c_str(), "%lf" , &psi[j][k]);
        if (i%NBasis==NBasis-1){k=-1;j++;}
    }
    k=0;
    file.clear();
    return 0;
}



int orb::read_orb_gamess( const char * nameorbs)
{
 ifstream in(nameorbs); // read in from the average file
 if (!in){
     cout << "Error file" << nameorbs << "  containg the orbitals non existant " <<endl;
     return 1;
 }
 string word;
 while( in ){
     in >> word;
     if ( word == "vectors" ){
        getline(in, word); //skip the rest of the line
        int i=0;
        int j=0;
        while ( i< NBasis){
            in >>word;
            sscanf(word.c_str(), "%lf", &psi[i][j]);
            j++;
            if(j==NBasis) {
                j=0;
                i++;
            }
        }
      }
   } 
 return 0;
}



void orb:: print_uhf_g03( const int & nel_A, const int & nel_B, const int & NBasis_A , const int NBasis_B) {
    //check nel even and NBasis+NBasis  == NBasis
    if ( nel_A % 2 != 0 || nel_B % 2 != 0 || NBasis_A + NBasis_B != NBasis ){
	cerr << "ERROR in input to print_uhf_g03" <<endl;
    }
	 
   
    FILE *out = fopen( "guess_card.g03", "w");	
    fprintf(out, "(1E15.8)\n");

    int count = 1;

    /* first the doubly occupied orbitals on A */
    for(int i = 0; i< nel_A/2-2 ; ++i) {
	fprintf(out, "\t%d Alpha", count);
	++count;
	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	

    /*then the double occupied orbitals on B */
    for (int i=NBasis_A; i < NBasis_A + nel_B/2-2 ; ++i){
	fprintf(out, "\t%d Alpha", count);
	++count;
    	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }
	
    /* then the singly occupied HOMO on A */
    fprintf(out, "\t%d Alpha", count);
    ++count;
    for (int j=0; j< NBasis ; j++ ){
	fprintf(out, "% 15.8E\n", psi[nel_A/2 -1][j]);
    }

    /* now the LUMOs of A */
    for(int i = nel_A; i< NBasis_A ; ++i) {
	fprintf(out, "\t%d Alpha", count);
	++count;
	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	
    /* and yes, the LUMOs of B*/
    for (int i=NBasis_A+nel_B; i < NBasis ; ++i){
	fprintf(out, "\t%d Alpha", count);
	++count;
    	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }

    fprintf(out,"\t0\n");

    /* and now the same thing but for the Beta electrons */
	
    /* first the doubly occupied orbitals on A */
    for(int i = 0; i< nel_A/2-2 ; ++i) {
	fprintf(out, "\t%d Beta", count);
	++count;
	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	

    /*then the double occupied orbitals on B */
    for (int i=NBasis_A; i < NBasis_A + nel_B/2-2 ; ++i){
	fprintf(out, "\t%d Beta", count);
	++count;
    	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }
	
    /* then the singly occupied HOMO on B */
    fprintf(out, "\t%d Beta", count);
    ++count;
    for (int j=0; j< NBasis ; j++ ){
	fprintf(out, "% 15.8E\n", psi[NBasis_A + nel_B - 1][j]);
    }

    /* now the LUMOs of A */
    for(int i = nel_A; i< NBasis_A ; ++i) {
	fprintf(out, "\t%d Beta", count);
	++count;
	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	
    /* and yes, the LUMOs of B*/
    for (int i=NBasis_A+nel_B; i < NBasis ; ++i){
	fprintf(out, "\t%d Beta", count);
	++count;
    	for (int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }


}



