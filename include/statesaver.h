/* 
 * File:   statesaver.h
 * Author: mayfalk
 *
 * Created on January 7, 2010, 3:29 PM
 */

#ifndef _STATESAVER_H
#define	_STATESAVER_H


#include <stdio.h>
#include "qmbead.h"
#include "qmtopology.h"
#include "qmpair.h"
#include "qmnblist.h"

using namespace std;

class StateSaver
{
public : 
    StateSaver(QMTopology & qmtop,string file);
    void Save();
    void Load();
    void Close();
    bool Seek(const int &);
      
    
    
private:
    void Write_PBCs();
    void Write_Molecules();
    void Write_QMBeads();
    void Write_QMNeighbourlist();

    void Read_PBCs();
    void Read_Molecules();
    void Read_QMBeads();
    void Read_QMNeighbourlist();

    template<typename T>
    void write(const T &v);

    template<typename T>
    T read(void);

    QMTopology *_qmtop;
    ofstream _out;
    ifstream _in;
    vector <streampos> _startpos;
};

template<typename T>
inline void StateSaver::write(const T &v)
{
  _out.write((const char *)&v, sizeof(T));
}

template<>
inline void StateSaver::write<std::string>(const string &v)
{
  write<unsigned short>(v.length());
  _out.write(v.c_str(), v.length());
}

template<typename T>
inline T StateSaver::read()
{
    T tmp;
    _in.read((char *)&tmp, sizeof(T));
    return tmp;
}

template<>
inline std::string StateSaver::read<std::string>()
{

  int L=read<unsigned short>();
  char tmp [L+1];
  _in.read(tmp, L);
  tmp[L]=0;
  return string(tmp);
}

#endif	/* _STATESAVER_H */
