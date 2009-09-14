/* 
 * File:   imcio.cpp
 * Author: ruehle
 *W
 * Created on September 14, 2009, 5:36 PM
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include "imcio.h"

typedef ub::symmetric_matrix<double> group_matrix;

void imcio_write_dS(string file, ub::vector<double> &r, ub::vector<double> &dS)
{
    // write the dS
    ofstream out_dS;
    out_dS.open(file.c_str());
    out_dS << setprecision(8);
    if (!out_dS)
        throw runtime_error(string("error, cannot open file ") + file);

    for (int i = 0; i < dS.size(); ++i) {
        out_dS << r[i] << " " << dS[i] << endl;
    }

    out_dS.close();
    cout << "written " << file << endl;
}

void imcio_write_matrix(string file, ub::symmetric_matrix<double> gmc)
{
    ofstream out_A;
    out_A.open(file.c_str());
    out_A << setprecision(8);

    if (!out_A)
        throw runtime_error(string("error, cannot open file ") + file);

    for (group_matrix::size_type i = 0; i < gmc.size1(); ++i) {
        for (group_matrix::size_type j = 0; j < gmc.size2(); ++j) {
            out_A << gmc(i, j) << " ";
        }
        out_A << endl;
    }
    out_A.close();
    cout << "written " << file << endl;
}

void imcio_write_index(string file, vector<string> &names, vector<int> &sizes)
{
    // write the index

    ofstream out_idx;
    out_idx.open(file.c_str());

    if (!out_idx)
        throw runtime_error(string("error, cannot open file ") + file);

    int last = 1;

    for (int i = 0; i < sizes.size(); ++i) {
        out_idx << names[i] << " " << last << ":" << last + sizes[i] - 1 << endl;
        last += sizes[i];
    }
    out_idx.close();
    cout << "written " << file << endl;
}
