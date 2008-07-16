#include <vector>
#include <string>
#include <iostream>

using namespace std;

//wrapper written with Denis to do MOO only on a file which give the unitary transformation + displacement one per line
void  parce_string (string line, string delims, vector<string>* result ) {
  string::size_type begIdx, endIdx;

  begIdx = line.find_first_not_of(delims);
  while (begIdx != string::npos) {
    endIdx = line.find_first_of (delims, begIdx);
    if (endIdx == string::npos) {
      endIdx = line.length();
    }
    result->push_back( line.substr(begIdx, endIdx-begIdx) );
    if (endIdx == line.length()) { break; cout << "I am still here";}
    begIdx = line.find_first_not_of (delims, endIdx);
  }
}


