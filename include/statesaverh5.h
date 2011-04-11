/* 
 * File:   statesaverh5.h
 * Author: ruehle
 *
 * Created on April 11, 2011, 2:00 PM
 */

#ifndef __VOTCA_MD2QM_STATESAVERH5_H
#define	__VOTCA_MD2QM_STATESAVERH5_H

#include <stdio.h>
#include "qmbead.h"
#include "qmtopology.h"
#include "qmpair.h"
#include "qmnblist.h"
#include "hdf5.h"
#include "hdf5_hl.h"

class StateSaverH5
{
public:
    StateSaverH5() {}
    ~StateSaverH5(){}

    void Open(QMTopology & qmtop,const string &file);
    void Close();
    void WriteFrame();
private:
    hid_t _file;
    int _frame;

    void WriteMolecules(hid_t context);
    void WriteBeads(hid_t context);
    void WritePairs(hid_t context);

    template<typename T>
    void write_attrib(hid_t location, string name, const T &v);

    template<typename T>
    void write_dataset(hid_t location, string name, const T &v);

    template<typename T>
    void write(hid_t location, string name, const T &v);


    QMTopology *_qmtop;
};

template<typename T>
void StateSaverH5::write(hid_t location, string name, const T &v)
{
    write_dataset<T>(location, name, v);
}


template<>
inline void StateSaverH5::write_attrib(hid_t location, string name, const double &v)
{
    H5LTset_attribute_double(location, ".", name.c_str(), &v, 1);
}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const double &v)
{
    hsize_t d=1;
    H5LTmake_dataset_double(location, name.c_str(), 1, &d, &v);

}

template<>
inline void StateSaverH5::write_attrib(hid_t location, string name, const int &v)
{
    H5LTset_attribute_int(location, ".", name.c_str(), &v, 1);
}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const int &v)
{
    hsize_t d=1;
    H5LTmake_dataset_int(location, name.c_str(), 1, &d, &v);

}

template<>
inline void StateSaverH5::write_attrib(hid_t location, string name, const string &v)
{
    H5LTset_attribute_string(location, ".", name.c_str(), v.c_str());
}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const string &v)
{
    H5LTmake_dataset_string(location, name.c_str(), v.c_str());
}

template<>
inline void StateSaverH5::write_attrib(hid_t location, string name, const vec &v)
{
    double tmp[3]={v.getX(), v.getY(), v.getZ()};
    H5LTset_attribute_double(location, ".", name.c_str(), tmp, 3);
}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const vec &v)
{
    double tmp[3]={v.getX(), v.getY(), v.getZ()};
    hsize_t d=3;
    H5LTmake_dataset_double(location, name.c_str(), 1, &d, tmp);
}

template<>
inline void StateSaverH5::write_attrib(hid_t location, string name, const matrix &v)
{
    matrix tmp=v;
    H5LTset_attribute_double(location, ".", name.c_str(), tmp[0], 9);
}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const matrix &v)
{
    hsize_t d=9;
    matrix tmp=v;
    H5LTmake_dataset_double(location, name.c_str(), 1, &d, tmp[0]);

}

template<>
inline void StateSaverH5::write_dataset(hid_t location, string name, const vector<double> &v)
{
    hsize_t d=v.size();
    H5LTmake_dataset_double(location, name.c_str(), 1, &d, &v[0]);
}


#endif	/* __VOTCA_MD2QM_STATESAVERH5_H */

