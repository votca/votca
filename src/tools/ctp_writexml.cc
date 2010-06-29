/* 
 * File:   playground.cpp
 * Author: mayfalk
 *
 * Created on May 21, 2010, 11:21 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include "write_xml.h"
#include "read_xml.h"

/*
 *
 */
int main(int argc, char** argv) {

    QMApplication qm_app;
    ReadXML reader;
    WriteXML writer;

    //qm_app.AddCalculator(dynamic_cast<QMCalculator*>(&reader));
    qm_app.AddCalculator(dynamic_cast<QMCalculator*>(&writer));
   
    qm_app.Run(argc, argv);

    return (EXIT_SUCCESS);
}

