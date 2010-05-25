/* 
 * File:   write_xml.h
 * Author: vehoff
 *
 * Created on May 25, 2010, 3:01 PM
 */

#ifndef _WRITE_XML_H
#define	_WRITE_XML_H

#include "qmcalculator.h"
#include "qmpair.h"

class WriteXML : public QMCalculator{
public:
    WriteXML() {};
    ~WriteXML() {};

    bool EvaluateFrame(int nr, int nframes, QMTopology *top);

private:
    ///  output stream to write integrals to xml file
    ofstream _out_int;
};

inline bool WriteXML::EvaluateFrame(int nr, int nframes, QMTopology *top){
    _out_int.open("integrals.xml", ios::app);
    QMNBList &nblist = top->nblist();
    if(_out_int!=0){
        if(nr==0){
            _out_int << "<qmtop>" << endl;
        }
        _out_int << "  <frame nr=\"" << nr << "\" />" << endl;
        for(QMNBList::iterator iter = nblist.begin();iter!=nblist.end();++iter){
            _out_int << "    <pair 1st=\"" << (*iter)->first->getId() << "\" 2nd=\"" << (*iter)->second->getId()
                    << "\" J_0=\"" << (*iter)->Js()[0]
                    << "\" Jeff=\"" << (*iter)->calcJeff2() << "\" />" << endl;
        }
        _out_int << "  </frame>" << endl;

        if (nr==(nframes-1)){
            _out_int << "</qmtop>" << endl;
        }
    }
    _out_int.close();
}

#endif	/* _WRITE_XML_H */

