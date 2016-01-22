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

#include <string>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

pair<double**, int> readS( string namefile)
{
    //read all the info from a punch file
    ifstream in(namefile.c_str());	 
    double **S;
    if (!in)
    {
	    cerr << "***Error*** : " << namefile << " doesn't exist.\n";
	    return make_pair(S,-1);
    }
    string word;	
    int Nbasis;
    while( in )
    {
	getline(in, word);
	if (   word == " Two-electron integral symmetry is turned off." 
	|| word == " Two-electron integral symmetry is turned on.")
	{
		in >> word;
		sscanf(word.c_str(), "%i", &Nbasis);
		cerr << "NBasis = " << Nbasis << endl;
	}

	if (word.find("*** Overlap")!=string::npos)
	{
	    S = new double * [Nbasis];

	    for (int m=0;m<Nbasis;m++)
	    {
		    S[m] = new double[Nbasis];		//get mem for over lap matrix
	    }	

	    for (int j=0; j<Nbasis; j+=5)			//go through the basis set in steps of 5
	    {
		    getline(in,word);				//Get labels 
		    int line=0;
		    for (int i=j; i<Nbasis; i++)
		    {
			    in >> word;        			//Get label
			    if (line<4)
			    {      			//Don't have all fields
				    for (int l=0; l<=line; l++)
				    {
					    in >> word;			
					    size_t pos=word.find("D");
					    word=word.substr(0,pos)+"e"+word.substr(pos+1);
					    sscanf(word.c_str(), "%lf" , &S[j+l][i]);
					    S[i][j+l]=S[j+l][i];
				    }
			    }
			    else
			    {					//Have 5 fields
				    for (int l=0; l<=4; l++)
				    {
					    in >> word;
					    size_t pos=word.find("D");
					    word=word.substr(0,pos)+"e"+word.substr(pos+1);
					    sscanf(word.c_str(), "%lf" , &S[j+l][i]);
					    S[i][j+l]=S[j+l][i];
				    }
			    }
		    line++;
		    }//for(i)
		    getline(in,word);				//Appears to be necessary to finish line
		}//for(j)
	    }//if(*** Overlap)

    }//while(in)
    cout << "Read S" << endl;
    return make_pair(S, Nbasis); 
}



int main(int argc, char **argv){
    if ( string(argv[1]).compare (  string("--help") ) == 0 ) {
	cout << "Usage: name log file, nameoutput " <<endl;	    
	return 1;
    } 

    pair <double**, int> res = readS(string(argv[1]));

    ofstream out(argv[2]);
    for (int i =0; i< res.second; i++){
    	for (int j=0; j< res.second; j++) {
	    out << scientific << setprecision(8) << res.first[i][j] << "  ";
	}
	out << '\n';
    }

    return 0;
}
