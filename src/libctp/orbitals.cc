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

#include "votca/ctp/orbitals.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace votca { namespace ctp {

 /* 
 * Reads in the MO coefficients from a GAUSSIAN fcheck file  
 */
bool Orbitals::ReadOrbitalsGaussian( const char * filename )
{
    
 ifstream _input_file ( filename ); 
 if ( _input_file.fail() ){
     cerr << "File "<< filename << " with molecular orbitals is not found " << endl;
     return 1;
 };

 string line;


 int i,j;
 int NBasis = 43;
 
 std::vector <string> file;
 std::vector <string> nrg;

 // number of coefficients per line is  given by the first line of the file (e.g. 5D15.8)
 int n_lin_in_wt = ( NBasis-1 ) / 5 + 2;
  
 int k = 0;
 i = 0;
 
// number of coefficients per line is  given by the first line of the file (e.g. 5D15.8)
 getline ( _input_file, line );
 std::vector<string> strs;
 boost::algorithm::split(strs, line, boost::is_any_of("(D)"));
 cout << strs.at(1) << endl;
 int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
 string format = strs.at(2);
 
 cout    << "Formatted orbital file: " 
         <<  nrecords_in_line << " records per line, in D" << format << " format." << endl;;

 unsigned _level = 0;
 
 while (_input_file){
	getline (_input_file, line);
        
        //cout << line << endl ;
        
        // if line has an equality sign, must be energy
        std::string::size_type energy_pos = line.find('=');
        
        if ( energy_pos != std::string::npos ) {
            cout << "energy line found at " << energy_pos << endl;            
            std::string energy_string = line.substr( energy_pos + 1 );
            cout << energy_pos << " " << energy_string << endl;
            /*energy_string.replace(energy_string.find("D"), sizeof("D")-1, "e");
            std::string::iterator _newend = std::remove(energy_string.begin(), energy_string.end(), ' ');
            energy_string.erase( _newend++, energy_string.end() );
            //double energy =  boost::lexical_cast<double>( energy_string );
            //cout << _level++ << " " << energy << endl; */
            
        } else {
                while (line.size() > 1 )
                {
                        /* string _coefficient;
                        _coefficient.assign(line,0,15);
                        cout << _coefficient << " " << endl;
                        
                        std::string::iterator _newend = std::remove(_coefficient.begin(), _coefficient.end(), ' ');
                        _coefficient.erase(_newend++,_coefficient.end());
                         cout << _coefficient << " ";
                        _coefficient.replace(_coefficient.find("D"), sizeof("D")-1, "e");
                        cout << _coefficient << " ";
                        //double coefficient = boost::lexical_cast<double>( _coefficient );
                        //cout << scientific << setprecision (7)<< coefficient << endl;
                        line.erase(0,15); */
                }
        }
   }
 
 
   if (file.size() != NBasis*NBasis ){
       throw runtime_error(string("I expect so many basis on this molecule ") + boost::lexical_cast<string>(NBasis*NBasis) +
               string(" but i read so many: " ) + boost::lexical_cast<string>(file.size()));
   }
   for(i=0,j=0,k=0;i<file.size();i++,k++)
   {

        //file[i].replace(file[i].size()-4,1,1,'e');
        //sscanf(file[i].c_str(), "%lf" , &psi[j][k]);
        //if (i%NBasis==NBasis-1){k=-1;j++;}
   }
   k=0;
   file.clear();
   return 0;
}


}}
