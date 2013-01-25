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


#include "idft.h"
#include <boost/algorithm/string.hpp>


namespace votca { namespace ctp {
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //
   
/* 
 * Reads in the MO coefficients from a GAUSSIAN fcheck file  
 */
bool IDFT::ReadOrbitalsGaussian( const char * filename )
{
    
 ifstream _input_file ( filename ); 
 if ( _input_file.fail() ){
     cerr << "File "<< filename << " with molecular orbitals is not found " << endl;
     return 1;
 };

 string line;
 string number;

 int i,j;
 int NBasis = 20;
 
 vector <string> file;
 vector <string> nrg;

 // number of coefficients per line is  given by the first line of the file (e.g. 5D15.8)
 int n_lin_in_wt = ( NBasis-1 ) / 5 + 2;
  
 int k = 0;
 i = 0;
 
// number of coefficients per line is  given by the first line of the file (e.g. 5D15.8)
 getline ( _input_file, line );
 vector<string> strs;
 boost::algorithm::split(strs, line, boost::is_any_of("(D)"));
 int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
 string format = strs.at(2);
 
 cout    << "Formatted orbital file: " 
         <<  nrecords_in_line << " records per line, D" << format << " format." << endl;;
 
 while (_input_file){
	getline (_input_file, line);
        cout << line << endl ;        
        
        if(k  % (n_lin_in_wt) != 0)
        {
                while (line.size() > 1 )
                {
                        number.assign(line,0,15);
                        cout << "DBG ONLY!!! " << number << endl;
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

// Line parser 
double IDFT::parsenrg(string & line){
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


/*
void IDFT::CleanUp() {

}
*/

/*
void IDFT::Initialize(Topology *top, Property *options) {

     cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;

    this->ParseOrbitalsXML(top, options);

    
}
*/

/*
void IDFT::ParseOrbitalsXML(Topology *top, Property *opt) {

}
*/

/*
void IDFT::EvalPair(Topology *top, QMPair *qmpair, int slot) {

}
*/

/*
void IDFT::CalculateJ(QMPair *pair) {

}
*/
    
    
}};