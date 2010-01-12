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
    StateSaver(QMTopology & qmtop);
    void Open(string file, bool bAppend = false);
    void Close();
    void Write_QMBeads(QMTopology *top);
    void Write_QMNeighbourlist(ofstream & out);
           
private:
    ofstream _out;
};

/*VICTOR:
class StateSaver
: public TrajectoryWriter
{
public:

    void Open(string file, bool bAppend = false);
    void Close();

    void RegisteredAt(ObjectFactory<string, TrajectoryWriter> &factory) {}

    void Write(Topology *conf);

private:
    FILE *_out;
};*/
#endif	/* _STATESAVER_H */
