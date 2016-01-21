/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/moo/orbitals.h>
#include <stdio.h>
#include <stdlib.h>

namespace votca { namespace moo {

int (orb::*orb::read_orb)(const char *)=&orb::read_orb_gauss;

void orb::init_orbitals (string * basis, const int & N, const char * namefile ){
    NBasis = N;
    bs.resize(NBasis);
    psi = new double* [NBasis];
    psi[0] = new double [NBasis * NBasis];
    bs[0]=basis[0];
    evl.resize(NBasis);
    for ( unsigned int i = 1 ; i < NBasis ; i ++){
            bs[i] = basis[i];
            psi[i] = psi[i-1] + NBasis;
    }
    (this->*read_orb) ( namefile);
}

void orb::reorder_for_libint(){
    /*reorder the d orbitals to agree with libint */
    unsigned int count=0;
    while (count < NBasis ){
    	if (bs[count] == "s") count +=1;
	else if (bs[count] == "x") count +=3;
	else if (bs[count] == "xx") { 
		bs[count+1] = "xy";
		bs[count+2] = "xz";
		bs[count+3] = "yy";	
		bs[count+4] = "yz";
		bs[count+5] = "zz";
		for (unsigned int j=0; j<NBasis; ++j){
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

void orb::strip_orbitals (const vector < int>& a){
    int nrorbs = a.size();
    double** psinew = new double*[nrorbs];
    psinew[0] = new double [nrorbs*NBasis];
    for(int i=1; i<nrorbs; i++){
        psinew[i] = psinew[i-1]+NBasis;
    }
    for (int i=0; i<nrorbs; i++){
        evl[i] = evl[a[i]];
        for (unsigned int j=0; j<NBasis; j++){
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
        for (unsigned int j=0; j<NBasis; j++){
            psi[i][j]=psinew[i][j];
        }
    }
    #ifdef DEBUG
    cout << "nrorbs: " << nrorbs << endl;
    cout << "psi[0][i] :" << endl;
    for (unsigned int j=0; j<NBasis; j++){
        cout << psi[0][j] << endl;
    }
    #endif
}

void orb::init_orbitals_stripped(const orb& orb1, const int& nrorbs ){
    NBasis = orb1.NBasis;
    bs.resize(NBasis);
    evl.resize(nrorbs);
    if(nrorbs==0) {
        psi = NULL;
    }
    else {
        psi = new double* [nrorbs];
        psi[0] = new double [nrorbs*NBasis];
        bs[0]=orb1.bs[0];
    }
    for(int i=1; i<nrorbs; i++){
        psi[i] = psi[i-1]+NBasis;
    }
    for(int i = 0; i < nrorbs; i++){
        evl[i] = orb1.evl[i];
        for ( unsigned int j = 0 ; j < NBasis ; j ++){
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

void orb::rot_orb(const unsigned int & orb , const double rot[3][3] ){
	
        matrix r(rot[0], rot[1], rot[2]);
        vector <pair <int, int> >::iterator itat=_basis_on_atom.begin();
        for ( ; itat!= _basis_on_atom.end(); ++itat){
            unsigned int first_basis = itat->first;
            unsigned int last_basis = first_basis + itat->second;
            unsigned int i=first_basis;
            while (i < last_basis){
             rot_orb(orb, &i, r);
            }
        }
}

inline void orb::rot_orb(const double rot[3][3]){
	for (unsigned int i =0; i< NBasis;i++){
		rot_orb( i, rot);
	}
}

void orb::rot_orb(const unsigned int & orb, unsigned int * i, const matrix &r){
    rot_orb(orb, i, psi[orb], r);
}

void orb::rot_orb(const unsigned int & orb , unsigned int *i, double * psi2, const matrix &r){
    //cout << "BS[" << *i << "] = " << bs[*i] << endl;
    //cout << "psi[orb][i] = psi[" << orb << "][" << *i <<  "]:" << psi[orb][*i] << endl;
    
    
            if ( bs[*i] == "s"){
                    (*i)++;
            }
            else if (bs[*i] == "x" ){
                    vec in (psi[orb][*i], psi[orb][*i+1],psi[orb][*i+2]) ;
                    vec out = r * in;
                    psi2[*i]=out.getX();
                    psi2[(*i)+1]=out.getY();
                    psi2[(*i)+2]=out.getZ();

                    /*double x_t,y_t,z_t;
                    psi2[*i] = psi[orb][*i]*rot.get(0,0) + psi[orb][*i+1] * rot.get(0,1) + psi[orb][*i+2] * rot.get(0,2);
                    psi2[*i+1]= psi[orb][*i]*rot.get(1,0) + psi[orb][*i+1] * rot.get(1,1) + psi[orb][*i+2] * rot.get(1,2);
                    psi2[*i+2]= psi[orb][*i]*rot.get(2,0) + psi[orb][*i+1] * rot.get(2,1) + psi[orb][*i+2] * rot.get(2,2) ;
                    */
                    (*i)+=3;
            }
            else if (bs[*i] =="xx" ){
                // fix a fixed reference frame
                    string internal_order[6] = {"xx", "xy", "xz", "yy", "yz", "zz"};
                    // find what order the basis sets are  actually in
                    vector <int> external_order;
                    for (int j =0; j<6;j++){
                        for(unsigned int k=*i;k<*i+6;k++){
                            if ( bs[k] == internal_order[j]){
                                external_order.push_back(k);
                                break;
                            }
                        }
                    }
                    double Dorbs[9] =
                    { psi[orb][external_order[0]],psi[orb][external_order[1]],psi[orb][external_order[2]],
                    psi[orb][external_order[1]],psi[orb][external_order[3]],psi[orb][external_order[4]],
                    psi[orb][external_order[2]],psi[orb][external_order[4]],psi[orb][external_order[5]]};
                    matrix in(Dorbs);

                    matrix out;
                    matrix tr(r);
		    matrix t(r);
                    out = t * in  * tr.Transpose();

                    psi2[external_order[0]] = out[0][0];
                    psi2[external_order[1]] = out[0][1];
                    psi2[external_order[2]] = out[0][2];
                    psi2[external_order[3]] = out[1][1];
                    psi2[external_order[4]] = out[1][2];
                    psi2[external_order[5]] = out[2][2];/*
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
                        ;*/

                    (*i)+=6;

            }
            else {
                    cerr << " Error in rot about axis" <<endl;
            }
}


void orb::rotate_someatoms(vector<int> atoms , matrix * M, 
        double * psi2 , const int & _orb){
        //cout << "orb = " << _orb << endl;
        vector <int>::iterator it;
        for (it = atoms.begin() ; it != atoms.end(); ++it){
            unsigned int first_basis = (_basis_on_atom[*it]).first;
            unsigned int last_basis = first_basis + (_basis_on_atom[*it]).second;
            unsigned int i=first_basis;
            while (i < last_basis){
                rot_orb(_orb, &i, psi2, *M); 
            }
        }
}

void orb::rot_orb(const unsigned int & orb , const matrix &rot){
	unsigned int i=0;
	while (i<NBasis)
	{
		rot_orb(orb, &i , rot);
	}
}

inline void orb::rot_orb(const matrix &rot){
	for (unsigned int i =0; i< NBasis;i++){
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

 unsigned int i,j;

 vector <string> file;
 vector <string> nrg;
 int n_lin_in_wt = (NBasis-1)/5+2;
 int k=0;
 i = 0;
 //skip the fortran specification line
 getline (in, line);
 while (in){
	getline (in, line);
        
        if(k  % (n_lin_in_wt) != 0)
        {
                while (line.size() > 1 )
                {
                        number.assign(line,0,15);
    //                    cout << "DBG ONLY!!! " << number << endl;
                        file.push_back(number);
                        line.erase(0,15);
                }
        }
        else {
            if (i <NBasis){
                evl[i]= parsenrg(line);
            }
            i++;
        }
        k++;
   }
   if (file.size() != NBasis*NBasis ){
       throw runtime_error(string("I expect so many basis on this molecule ") + boost::lexical_cast<string>(NBasis*NBasis) +
               string(" but i read so many: " ) + boost::lexical_cast<string>(file.size()));
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
        unsigned int i=0;
        unsigned int j=0;
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


void orb::print_g03(string & name, string  mode){
	FILE *out = fopen( name.c_str(), mode.c_str());	
	fprintf(out, "(5E15.8)\n");
        int count=1;
	for (unsigned int i=0;i < NBasis;i++){
		fprintf(out, "\t%d Alpha\n", count);
		++count;
		for (unsigned int j=0; j< NBasis ; j++ ){
			fprintf(out, "% 15.8E", psi[i][j]);
                        if ((j+1)%5==0 || j==NBasis-1) fprintf(out,"\n");
		}
	}
        fprintf(out, "\n");
        fclose(out);
}

void orb:: print_uhf_g03( const int & nel_A, const int & nel_B, const unsigned int & NBasis_A , const unsigned int NBasis_B) {
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
	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	

    /*then the double occupied orbitals on B */
    for (unsigned int i=NBasis_A; i < NBasis_A + nel_B/2-2 ; ++i){
	fprintf(out, "\t%d Alpha", count);
	++count;
    	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }
	
    /* then the singly occupied HOMO on A */
    fprintf(out, "\t%d Alpha", count);
    ++count;
    for (unsigned int j=0; j< NBasis ; j++ ){
	fprintf(out, "% 15.8E\n", psi[nel_A/2 -1][j]);
    }

    /* now the LUMOs of A */
    for(unsigned int i = nel_A; i< NBasis_A ; ++i) {
	fprintf(out, "\t%d Alpha", count);
	++count;
	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	
    /* and yes, the LUMOs of B*/
    for (unsigned int i=NBasis_A+nel_B; i < NBasis ; ++i){
	fprintf(out, "\t%d Alpha", count);
	++count;
    	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }

    fprintf(out,"\t0\n");

    /* and now the same thing but for the Beta electrons */
	
    /* first the doubly occupied orbitals on A */
    for(int i = 0; i< nel_A/2-2 ; ++i) {
	fprintf(out, "\t%d Beta", count);
	++count;
	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	

    /*then the double occupied orbitals on B */
    for (unsigned int i=NBasis_A; i < NBasis_A + nel_B/2-2 ; ++i){
	fprintf(out, "\t%d Beta", count);
	++count;
    	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }
	
    /* then the singly occupied HOMO on B */
    fprintf(out, "\t%d Beta", count);
    ++count;
    for (unsigned int j=0; j< NBasis ; j++ ){
	fprintf(out, "% 15.8E\n", psi[NBasis_A + nel_B - 1][j]);
    }

    /* now the LUMOs of A */
    for(unsigned int i = nel_A; i< NBasis_A ; ++i) {
	fprintf(out, "\t%d Beta", count);
	++count;
	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }	
    /* and yes, the LUMOs of B*/
    for (unsigned int i=NBasis_A+nel_B; i < NBasis ; ++i){
	fprintf(out, "\t%d Beta", count);
	++count;
    	for (unsigned int j=0; j< NBasis ; j++ ){
		fprintf(out, "% 15.8E\n", psi[i][j]);
	}
    }


}




void orb::dimerise_orbs(const orb & A, const orb & B, const int &elA, const int &elB) {
    if ( psi != 0) {
        clear();
    }
    NBasis = A.NBasis + B.NBasis;

    /* set up the arrays */
    bs.resize(NBasis);
    evl.resize(NBasis);

    if(NBasis==0)
        psi=NULL;
    else {
        psi = new double* [NBasis];
        psi[0] = new double [NBasis * NBasis];
    }
    for ( unsigned int i = 1 ; i < NBasis ; i ++){
            psi[i] = psi[i-1] + NBasis;
    }

    /* cp bs info */
    for (unsigned int i=0; i< A.NBasis ;++i ){
        bs[i] = A.bs[i];
        evl[i] =0.;
    }
    for (unsigned int i=0; i< B.NBasis ;++i ){
        bs[A.NBasis + i ] = B.bs[i];
        evl[A.NBasis + i ]=0.;
    }

    /* copy occupied orbitals*/
    int occA = (elA)/2;
    int occB = (elB)/2;
    
    /*copy orbitals  from A*/
    for (int i=0; i <occA;++i){
        for (unsigned int j=0 ; j< A.NBasis ;++j){
            psi[i][j] = A.psi[i][j];
        }
        for (unsigned int j=0; j< B.NBasis ;++j){
            psi[i][A.NBasis+j] =0.0;
        }
    }

    /*copy orbitals from B*/
    for (int i=0; i <occB ;++i){
        for (unsigned int j=0; j< A.NBasis ;++j){
            psi[occA+i][j] = 0.0;
        }
        for (unsigned int j=0; j< B.NBasis ;++j){
            psi[occA+i][A.NBasis+j] = B.psi[i][j];
        }
    }

    /*now copy the empty ones*/

    /*copy orbitals  from A*/
    for (unsigned int i=occA; i <A.NBasis;++i){
        for (unsigned int j=0 ; j< A.NBasis ;++j){
            psi[occB+i][j] = A.psi[i][j];
        }
        for (unsigned int j=0; j< B.NBasis ;++j){
            psi[occB+i][A.NBasis+j] =0.0;
        }
    }

    /*copy orbitals from B*/
    for (unsigned int i=occB; i <B.NBasis ;++i){
        for (unsigned int j=0; j< A.NBasis ;++j){
            psi[A.NBasis+i][j] = 0.0;
        }
        for (unsigned int j=0; j< B.NBasis ;++j){
            psi[A.NBasis+i][A.NBasis+j] = B.psi[i][j];
        }
    }

}

double orb::parsenrg(string & line){
    size_t find = line.find_last_of("=");
    if (find >= line.length() ){
        return 0.;
    }
    string found = line.substr(find+1, line.length() - find);
    double r;
    found.replace(found.size()-4,1,1,'e');
    sscanf(found.c_str(), "%lf" , &r);

    return r;
}

}}
